#' Take a RangedSummarizedExperiment x & turn it into a MultiAssayExperiment y
#' where the former metadata(x) becomes the integrity-checked linkedAssays(y).
#'
#' @param x   a RangedSummarizedExperiment whose metadata() contains assays
#' 
#' @return    a MultiAssayExperiment y with integrity-checked linkedAssays(y)
#' 
#' @export
#' 
convertToMultiAssay <- function(x) { 
  stopifnot(is(x, "RangedSummarizedExperiment"))
  x <- as(x, "MultiAssayExperiment")
  x@perSampleMetadata <- lapply(x@metadata, function(y) y[, colnames(x)])
  x@metadata <- list()
  return(x)
}
