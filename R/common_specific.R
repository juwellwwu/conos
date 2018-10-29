
common_specific <- function(conObj, groups=NULL, cooksCutoff = FALSE, min.cell.count = 10,
                            independentFiltering = FALSE, n.cores=1, cluster.sep.chr = '+',return.details=TRUE,
                            exclude.samples = c(), sign.cutoff = 0.05, use.corrected.pval = FALSE, 
                            ct.specific.thres = 1, ct.common.thres = 5 ) {
  
  ## DEVEL
  min.cell.count <- 10
  conObj <- con
  groups <- as.factor(con$clusters$`multi level`$groups)
  cluster.sep.chr <- '+'
  independentFiltering = FALSE
  cooksCutoff = FALSE
  n.cores <- 4
  exclude.samples <- c()
  ###
  sign.cutoff <- 0.05
  use.corrected.pval <- FALSE
  ct.specific.thres <- 1
  ct.common.thres <- 5
  ##
  
  ## Use all but excluded samples
  samples.used <- names(conObj$samples)[!names(conObj$samples) %in% exclude.samples]
  
  ## Generate an aggregated matrix
  raw.mats <- lapply(conObj$samples[samples.used], function(p2) {
    p2$misc$rawCounts
  })
  
  common.genes <- Reduce(intersect,lapply(raw.mats, colnames))
  raw.mats <- lapply(raw.mats, function(x) {x[,common.genes]})
  aggr2 <- lapply(raw.mats, function(x) {
    g1 <- groups[intersect(names(groups), rownames(x))]
    t1 <- as.numeric(table(g1))
    names(t1) <- levels(g1);
    droplevels <- names(t1)[t1 < min.cell.count]
    g1.n <- names(g1)
    g1 <- as.character(g1)
    names(g1) <- g1.n
    g1[g1 %in% droplevels] <- NA
    g1 <- as.factor(g1)
    aggr <- Matrix.utils::aggregate.Matrix(x, g1)
    aggr <- aggr[rownames(aggr) != "NA",]
    aggr
  })
  
  aggr2 <- lapply(names(aggr2), function(n) {
    x <- aggr2[[n]]
    rownames(x) <- paste0(n,cluster.sep.chr,rownames(aggr2[[n]]))
    x
  })
  
  aggr2 <- t(do.call(rbind, aggr2))
  rm(raw.mats); gc()
  
  aggr2.meta <- data.frame(
    row.names = colnames(aggr2),
    sample = colnames(aggr2),
    library = factor(strpart(colnames(aggr2),'+',1,fixed='T')),
    celltype = factor(as.character(strpart(colnames(aggr2),'+',2,fixed='T')))
  )
  
  ## For every individual perform differential expression for each celltype between the current sample and everything else
  de.res <- parallel::mclapply(nbHelpers::namedLevels(aggr2.meta$library), function(lib) {
    ## For each individual
    de.res.1 <- lapply(nbHelpers::namedLevels(aggr2.meta$celltype), function(ct) {
      ## For each cell type
      try({
        ##dev
        #ct <- nbHelpers::namedLevels(aggr2.meta$celltype)[1]
        #lib <- nbHelpers::namedLevels(aggr2.meta$library)[1]
        tmp.meta <- aggr2.meta[aggr2.meta$celltype == ct,]
        tmp.meta$compVector <- factor(ifelse(tmp.meta$library == lib,'cursample','other'))
        tmp.meta$compVector <- relevel(tmp.meta$compVector, ref = 'other')
        dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = aggr2[,rownames(tmp.meta)], colData = tmp.meta, design = ~ compVector)
        dds1 <- DESeq2::DESeq(dds1)
        res1 <- DESeq2::results(dds1, cooksCutoff = cooksCutoff, independentFiltering = independentFiltering)
        res1 <- as.data.frame(res1)
        res1 <- res1[order(res1$padj,decreasing = FALSE),]
        res1
      })
    })
    de.res.1
  }, mc.cores=n.cores)
  
  ## config
  # sign.cutoff <- 0.05
  # use.corrected.pval <- FALSE
  # ct.specific.thres <- 1
  # ct.common.thres <- 5
  
  ## What is common between the cell different cell types for each individual?
  per.individual <- parallel::mclapply(de.res, function(lib) {
    all.ct.genes <- Reduce(union,lapply(lib, rownames))
    per.ct.genes <- do.call(cbind, lapply(lib, function(x) {
      if(use.corrected.pval) {
        score.val <- x$padj
      } else {
        score.val <- x$pvalue
      }
      all.ct.genes %in% rownames(x)[score.val < sign.cutoff]
    }))
    colnames(per.ct.genes) <- names(lib)
    rownames(per.ct.genes) <- all.ct.genes
    no.times.gene.sign <- apply(per.ct.genes, 1, sum)
    specific.genes <- names(no.times.gene.sign)[no.times.gene.sign <= ct.specific.thres]
    common.genes <- names(no.times.gene.sign)[no.times.gene.sign >= ct.common.thres]
    list(specific.genes = specific.genes, common.genes = common.genes)
  }, mc.cores = n.cores)
  
  ## Per cluster specific genes
  
  
  ## Return
  list(per.individual = per.individual)
}