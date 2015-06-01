## semi-automated DMRcation with backoff 
## FIXME: allow automating survival DMRs 
## 
processDMRs <- function(grSet, name, design=NULL, ...) {

  dmrFile <- paste0(name, "_DMRs.rds")
  if (!is.null(design)) {
    message("Calling DMRs from design matrix...")
  } else if (!any(grSet$Sample_Group == "control")) {
    message("Calling DMRs against controls...")
    design <- model.matrix( ~ as.numeric(grSet$Sample_Group != "control"))
  } else { 
    message("No design matrix and no controls provided.  Skipping DMR calling.")
    invisible(NULL)
  }
  DMRs <- getDMRs(grSet, design=design, ...)
  saveRDS(DMRs, dmrFile)
  DMRsGR <- extractRanges(DMRs) 
  export(DMRsGR, paste0(name, "_DMRs.hg19.bed"))
  message("Done.")
  invisible(DMRsGR)

} 

