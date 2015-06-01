## components of a simple 450k + CNV + expression data pipeline: CopyNumber450k
##
processCNV <- function(rgSet, name, ...) {
  ## process the CNV aspect of the data 
  cnvSetFile <- paste0(name, "_CNV.rds")
  samples <- sampleNames(rgSet)
  if (!any(rgSet$Sample_Group == "control")) { ## EISAI, Klco, Birmingham
    library(CopyNumber450kData)
    data(RGcontrolSetEx)
    rgSet <- combine(rgSet, RGcontrolSetEx)
  }
  cnvSet <- CNV450kSet(rgSet)
  cnvSet <- computeSignificance(segmentize(normalize(dropSNPprobes(cnvSet))))
  cnvSet <- cnvSet[ , samples] ## get rid of excess controls 
  message("Saving CNV results to ", cnvSetFile, "...")
  saveRDS(cnvSet, file=cnvSetFile)
  message("Saving segmentized CNV calls to ", 
          sub("rds", "signif.hg19.csv", cnvSetFile), "...")
  write.csv(cnvSet, file=sub("rds", "signif.hg19.csv", cnvSetFile))
  postprocessCNV(name)
}
