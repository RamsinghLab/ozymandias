## funnorm + pvals + betas + SNPs 
processMeth <- function(rgSet, name, pcutoff=0.01, saveExcel=FALSE, ...) {

  grSet <- preprocessFunnorm(rgSet)
  metadata(grSet)$SNPs <- getSnpBeta(rgSet)
  assays(grSet)$pval <- matrix(NA_real_, ncol=ncol(grSet), nrow=nrow(grSet))
  assays(grSet)$pval <- detectionP(rgSet)[rownames(grSet), colnames(grSet)]
  is.na(assays(grSet)$Beta[which(assays(grSet)$pval > pcutoff)]) <- TRUE 

  if (saveExcel == TRUE) {
    betaFile <- paste0(name, "_betas.xls")
    message("Saving processed beta values to ", betaFile, "...")
    write.csv(getBeta(grSet), betaFile)

    pvalFile <- paste0(name, "_pvals.xls")
    message("Saving processed pvalues to ", pvalFile, "...")
    write.csv(assays(grSet)$pval, pvalFile)
  }

  grSetFile <- paste0(name, "_450k.rds")
  message("Saving GenomicRatioSet to ", grSetFile, "...")
  saveRDS(grSet, file=grSetFile)
  invisible(grSet)

}
