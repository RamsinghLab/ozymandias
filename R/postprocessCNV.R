postprocessCNV <- function(name) { 

  ## take the results from CopyNumber450k and process them more usefully
  message("Reformatting significant segments into GRanges/BED format...")
  sigCNV <- read.csv(paste0(name, "_CNV.signif.hg19.csv"),
                     row.names=1, stringsAsFactors=FALSE)
  sigCNV <- sigCNV[which(sigCNV$isSignificant == TRUE), ]
  names(sigCNV) <- sub("loc.start", "chromStart", names(sigCNV))
  names(sigCNV) <- sub("loc.end", "chromEnd", names(sigCNV))
  sigCNV$score <- sigCNV$seg.mean
  sigCNV$name <- paste0(sigCNV$Sample, ":", 
                        "adj.p=", round(sigCNV$adjusted.pvalue, 6))
  sigCNV$what <- ifelse(sigCNV$seg.mean > 0, "gain", "loss")
  sigCNV <- sigCNV[, c("chrom","chromStart","chromEnd","score","name","what")]

  message("Constructing a GRanges from significant gains & losses...")
  data(seqinfo.hg19) ## since we have to do this on the regular...
  cnvGR <- makeGRangesFromDataFrame(sigCNV, keep.extra.columns=TRUE)
  seqinfo(cnvGR) <- seqinfo.hg19[seqlevels(cnvGR)]
  saveRDS(cnvGR, file=paste0(name, "_CNV.GRanges.hg19.rds"))

  message("Dumping significant CNV segments to a BED file...")
  export(cnvGR, paste0(name, "_CNV.signif.hg19.bed"))
  
  message("Dumping significant copy number losses to a BED file...")
  export(split(cnvGR, cnvGR$what)$loss, paste0(name, "_CN_LOSS.hg19.bed"))

  message("Dumping significant copy number gains to a BED file...")
  export(split(cnvGR, cnvGR$what)$gain, paste0(name, "_CN_GAIN.hg19.bed"))

  ## to tally recurrences: 
  getDepth <- function(x, what) { 
    stopifnot(what %in% c("loss","gain"))
    bins <- disjoin(x)
    bins$score <- countOverlaps(bins, x)
    if (what == "loss") bins$score <- -1 * bins$score
    return(bins)
  }

  ## depth of recurrence: losses
  message("Dumping focal depth of copy number losses to a BigWig file...")
  export(getDepth(split(cnvGR, cnvGR$what)$loss, "loss"),
         paste0(name, "_CN_LOSS_DEPTH.hg19.bw"))

  ## depth of recurrence: gains
  message("Dumping focal depth of copy number gains to a BigWig file...")
  export(getDepth(split(cnvGR, cnvGR$what)$gain, "gain"),
         paste0(name, "_CN_GAIN_DEPTH.hg19.bw"))

  invisible(cnvGR)

}
