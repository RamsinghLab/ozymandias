#' Dump DNA methylation values, gene-level counts, or whatever for IGV viewing
#'
#' @param x         a SummarizedExperiment-like object of some sort
#' @param assay     the name of the assay to dump (default: assays(x)[[1]])
#' @param group     a grouping factor to dump by (default: dump everyone) 
#' @param stub      the filename stub to dump (default: name of x, or x_group)
#' @param genome    the genome to use (default is hg19, due to annotations)
#' 
#' @return list     a list of filenames dumped
#'
asTDF <- function(x, assay=1, group=NULL, stub=NULL, genome="hg19") { 
  if (is.null(stub)) stub <- as.character(match.call()["x"])
  asGCT(x, toTDF=T, assay=assay, group=group, stub=stub, genome=genome)
}

.toTDF <- function(filename, probeGR) { # {{{
  fileout <- sub("gct", "tdf", filename)
  probefile <- paste0(sub(".gct", "", filename), ".probes.bed")
  .toBED4(probeGR, probefile)
  command <- paste("igvtools", "totdf", 
                   "-p", probefile, 
                   filename,
                   fileout,
                   "hg19")
  message("Converting to TDF via:")
  message(command)
  retval <- system(command)
  if (file.exists(fileout)) {
    unlink(probefile)
    unlink(filename)
    return(fileout)
  } else {
    return(filename)
  }
} # }}}
