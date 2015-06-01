## for tidying up rowRanges prior to merging
##
stripMcols <- function(x) { # {{{
  if (is(x, "GRangesList")) {
    xx <- unlist(x)
    nm <- unlist(lapply(strsplit(names(xx), ".", fixed=T), `[`, 1))
    mcols(xx) <- NULL
    x <- split(xx, nm)
  } else if (is(x, "GRanges")) {
    mcols(x) <- NULL
  }
  return(x)
} # }}}

## mask beta values by pval for a GenomicRatioSet 
## 
maskBetas <- function(grSet, pcutoff=0.01, ...) { # {{{

  if (!"pval" %in% names(assays(grSet))) {
    message("You do not have detection calls stored in assays(grSet)$pval")
    message("Cannot mask beta values.  Returning your grSet unaltered.")
  } else if (!is(assays(grSet)$pval, "Matrix")) { 
    ## this can save a HUGE amount of space overall 
    assays(grSet)$pval <- Matrix(assays(grSet)$pval)
  }
  is.na(assays(grSet)$Beta[which(assays(grSet)$pval > pcutoff)]) <- TRUE 
  return(grSet)

} # }}}
