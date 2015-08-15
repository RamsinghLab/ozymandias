#' Pull a 27k dataset from GEO.
#' 
#' Doesn't matter what the annotations are or any other bullshit. Just get it.
#' 
#' @param GSE   The dataset.
#'
#' @return a GenomicRatioSet.
#'
fetch27k <- function(GSE) {

  gset <- getGEO(GSE)[[1]]
  library(FDb.InfiniumMethylation.hg19)
  hm27 <- get27k()
  GenomicRatioSet(gr=hm27[featureNames(gset)], 
                  pData=pData(gset),
                  Beta=exprs(gset)) 
}
