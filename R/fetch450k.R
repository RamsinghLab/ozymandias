#' Pull a 450k dataset from GEO.
#' 
#' Doesn't matter what the annotations are or any other bullshit. Just get it.
#' 
#' @param GSE   The dataset.
#'
#' @return a GenomicRatioSet.
#'
fetch450k <- function(GSE) {
  gset <- getGEO(GSE)[[1]]
  library(FDb.InfiniumMethylation.hg19)
  hm450 <- get450k()
  prepForReanalysis(GenomicRatioSet(gr=hm450[featureNames(gset)], 
                                    pData=pData(gset),
                                    Beta=exprs(gset)))
}
