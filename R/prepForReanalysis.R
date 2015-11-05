#' extract IDAT names, etc. for reprocessing from raw data
#' 
#' @param res   a SummarizedExperiment, GenomicRatioSet, or suchlike
#' 
#' @return      the same thing but with tidied-up Basename, Sample_Name, etc.
#' 
#' @export
#'
prepForReanalysis <- function(res) {
  if (!is.null(res$supplementary_file)) {
    res$Sample_Name <- res$title
    res$Sample_Group <- "Experimental"
    res$Basename <- gsub(".gz","", 
                         gsub("_Grn.idat","", 
                              basename(res$supplementary_file)))
  }
  return(res)
}
