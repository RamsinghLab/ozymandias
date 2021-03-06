#' Dump DNA methylation values, gene-level counts, or whatever for IGV viewing
#'
#' @param x         a SummarizedExperiment-like object of some sort
#' @param toTDF     run igvtools totdf on the .GCT files? (default: no)
#' @param assay     the name of the assay to dump (default: assays(x)[[1]])
#' @param group     a grouping factor to dump by (default: dump everyone) 
#' @param stub      the filename stub to dump (default: name of x, or x_group)
#' @param genome    the genome to use (default is hg19, due to annotations)
#' 
#' @return list     a list of filenames dumped
#'
asGCT <- function(x, toTDF=F, assay=1, group=NULL, stub=NULL, genome="hg19") { 

  x <- sort(x) ## must be sorted by position for GCT and toTDF 
  if (is.null(stub)) stub <- as.character(match.call()["x"])
  if (is.null(assay)) assay <- 1 

  if (is.null(group)) { 
    ## dump everyone at the same time 
    filename <- .dumpGCT(assays(x)[[assay]], stub, genome)
    if (toTDF) filename <- .toTDF(filename, SummarizedExperiment::rowRanges(x))
    filenames <- list(stub=filename)
  } else { 
    ## dump by group
    filenames <- list()
    for (g in levels(as.factor(group))) {
      grp <- which(group == g)
      grpstub <- paste(stub, g, sep=".")
      filename <- .dumpGCT(assays(x[, grp])[[assay]], grpstub, genome)
      if (toTDF) {
        filename <- .toTDF(filename, probeGR=SummarizedExperiment::rowRanges(x))
      }
      filenames <- append(filenames, filename)
    }
  }

  invisible(filenames)

}

.dumpGCT <- function(mat, stub, genome="hg19") { # {{{
  fname <- paste(stub, genome, "gct", sep=".")
  contents <- data.frame(NAME=rownames(mat), Description=rownames(mat), mat)
  message("Dumping values to ", fname, "...")
  ## formatting for GCT 
  cat("#1.2", "\n", file=fname)
  cat(nrow(mat), "\t", ncol(mat), "\n", file=fname, append=TRUE) 
  cat(colnames(mat), "\n", file=fname, sep="\t", append=TRUE) 
  write.table(contents,
              file=fname, append=T, row.names=F, col.names=F, quote=F, sep="\t")
  return(fname)
} # }}}

.toBED4 <- function(gr, filename) { # {{{
  coord <- as.data.frame(gr)[,1:3]
  coord[,2] <- coord[,2] - 1 ## derrrrp
  coord <- cbind(coord, rownames(coord))
  write.table(coord, row.names=F, col.names=F, quote=F, file=filename, sep="\t")
} # }}}
