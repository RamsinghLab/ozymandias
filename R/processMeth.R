## funnorm + pvals + betas + SNPs 
processMeth <- function(rgSet, name, pcutoff=0.01, ...) {

  grSet <- preprocessFunnorm(rgSet)
  assays(grSet)$pval <- detectionP(rgSet)
  metadata(grSet)$SNPs <- getSnpBeta(rgSet)
  is.na(assays(grSet)$Beta[which(assays(grSet)$pval > pcutoff)]) <- TRUE 

  betaFile <- paste0(name, "_betas.xls")
  message("Saving processed beta values to ", betaFile, "...")
  write.csv(getBeta(grSet), betaFile)

  pvalFile <- paste0(name, "_pvals.xls")
  message("Saving processed pvalues to ", pvalFile, "...")
  write.csv(assays(grSet)$pval, pvalFile)

  grSetFile <- paste0(name, "_450k.rds")
  message("Saving GenomicRatioSet to ", grSetFile, "...")
  saveRDS(grSet, file=grSetFile)
  invisible(grSet)

}
