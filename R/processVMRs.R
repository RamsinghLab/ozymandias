## this could be extended to incorporate WGBS and/or TFBS scanning if need be 
##
processVMRs <- function(grSet, name, ...) {

  message("Calling VMRs...")
  vmrFile <- paste0(name, "_VMRs.rds")
  VMRs <- getVMRs(grSet, ...)
  saveRDS(VMRs, vmrFile)
  VMRsGR <- extractRanges(VMRs)
  export(VMRsGR, paste0(name, "_VMRs.hg19.bed"))
  message("Done.")
  invisible(VMRsGR)
  
}
