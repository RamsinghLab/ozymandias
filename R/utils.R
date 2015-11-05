#' array-ify a string with a possible leading path/URL/whatever
#'
#' @param x   a string
#' @param y   a split character
#' 
#' @return    the array produced by splitting basename(x) on y
#'
#' @export
#'
str2vec <- function(x, y="_") {
  strsplit(basename(x), y)[[1]]
}

#' grab the z'th element(s) off of a vector formed by splitting x on y 
#' useful for (e.g.) grabbing source barcodes from GEO 
#' 
#' @param x   a string
#' @param y   a split character ("_")
#' @param z   which element(s) to retrieve (1)
#' 
#' @return    the element(s) [z] produced by splitting basename(x) on y
#'
#' @export
#'
strshift <- function(x, y="_", z=1) {
  paste(str2vec(x, y)[z], collapse=y)
}

#' pop a string off a string, for (e.g.) grabbing source barcodes from GEO 
#' 
#' @param x   a string
#' @param y   a split character
#' @param z   which element(s) to retrieve (the last one)
#' 
#' @return    the element(S) produced by splitting basename(x) on y
#'
#' @export
#'
strpop <- function(x, y="_", z=NULL) {
  vec <- str2vec(x, y)
  if (is.null(z)) z <- length(vec)
  paste(vec[z], collapse=y)
}

#' The opposite of strpop: push string(s) onto a string.
#' 
#' @param w   a string
#' @param x   one or more elements to push onto the string
#' @param y   optional, split character (default: "_") 
#' @param z   optional, add x after str2vec(w)[z] (default: length(str2vec(w)))
#' 
#' @return    a string
#'
#' @export
#'
strpush <- function(w, x, y="_", z=NULL) {
  vec <- str2vec(w, y)
  if (is.null(z)) z <- length(vec)
  paste(append(vec, x, z), collapse=y)
}

#' for tidying up rowRanges prior to merging
#'
#' @param x  a GRanges or GRangesList that has metadata columns it doesn't need 
#'
#' @describeIn  ozymandias-utils
#'
#' @export
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

#' mask (set to NA) beta values by pval for a GenomicRatioSet 
#' 
#' @param grSet   a GenomicRatioSet with matrices "Beta" and "pval"
#' @param pcutoff optional p-value cutoff for masking, default is 0.01
#' 
#' @describeIn  ozymandias-utils
#'
#' @export
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

#' convenience function to feed getSurv, OSurv, EFSurv
#'
#' @param x         a SummarizedExperiment-like dataset
#' @param timeCol   column name for time to event covariate, default "time"
#' @param eventCol  column name for censoring/event covariate, default "event"
#'
#' @return          a boolean vector: whether each column has both covariates
#'
#' @describeIn  ozymandias-utils
#'
#' @export
hasSurv <- function(x, timeCol="time", eventCol="event") {
  return(!is.na(x[[timeCol]]) & !is.na(x[[eventCol]]))
}

#' convenience function to feed OSurv, EFSurv
#'
#' @param x         a SummarizedExperiment-like dataset
#' @param timeCol   column name for time to event, default "OS"
#' @param eventCol  column name for censoring/event (0/1), default "OSevent"
#'
#' @return          a \code{Surv} object, if one can be constructed from x 
#'
#' @describeIn  ozymandias-utils
#'
#' @export
getSurv <- function(x, timeCol="time", eventCol="event") { 
  xx <- colData(x)[hasSurv(x, timeCol=timeCol, eventCol=eventCol),]
  Surv(xx[[timeCol]], xx[[eventCol]])
}

#' Get overall survival for a SummarizedExperiment-like dataset
#' 
#' @param x         a SummarizedExperiment-like dataset
#' @param timeCol   column name for time to event, default "OS"
#' @param eventCol  column name for censoring/event (0/1), default "OSevent"
#' 
#' @return          a \code{Surv} object, if one can be constructed from x 
#'
#' @describeIn  ozymandias-utils
#'
#' @export
OSurv <- function(x, timeCol="OS", eventCol="OSevent") {
  getSurv(x, timeCol, eventCol)
}

#' Get event-free survival from SummarizedExperiment-like data
#' 
#' @param x         a SummarizedExperiment-like dataset
#' @param timeCol   column name for time to event, default "EFS"
#' @param eventCol  column name for censoring/event (0/1), default "EFSevent"
#' 
#' @return          a \code{Surv} object, if one can be constructed from x 
#'
#' @describeIn  ozymandias-utils
#'
#' @export
EFSurv <- function(x, timeCol="EFS", eventCol="EFSevent") {
  getSurv(x, timeCol, eventCol)
}
