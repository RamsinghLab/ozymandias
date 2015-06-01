## this is kind of silly but it encapsulates two important defaults.
##
processIDATs <- function(targets, name, ...) {
  rgSetFile <- paste0(name, "_rgSet.rds")
  ## failures can be mysterious if we do not set verbose!
  rgSet <- read.450k.exp(base=".", targets, verbose=TRUE)
  sampleNames(rgSet) <- targets$Sample_Name ## mandatory!
  message(ncol(rgSet), " samples processed sucessfully.")
  message("Saving raw intensities to ", rgSetFile, "...")
  saveRDS(rgSet, rgSetFile)
  invisible(rgSet)
}

