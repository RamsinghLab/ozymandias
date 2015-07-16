#' Something like a SummarizedExperiment, but with "holes" and metadata checks
#'
#' @slot perSampleMetadata associated per-sample data, checked against colnames
#'
#' @export
#'
setClass("MultiAssayExperiment", 
         representation(perSampleMetadata="list"), # or a DataFrameList maybe?
         contains="RangedSummarizedExperiment") # RangedSummarizedExperiments?

setGeneric("perSampleMetadata", 
           function(object) standardGeneric("perSampleMetadata"))

setMethod("perSampleMetadata", "MultiAssayExperiment", 
           function(object) object@perSampleMetadata)

setMethod("show", "MultiAssayExperiment", 
          function(object) {
            callNextMethod()
            some <- function(x) { 
              if(length(x) %in% 1:3) {
                return(x) 
              } else {
                c(x[1], "...", x[length(x)])
              }
            }
            cat("perSampleMetadata:", names(object@perSampleMetadata), "\n")
            for (i in names(perSampleMetadata(object))) { 
              cat(paste0("  colnames(perSampleMetadata(", i, ")): "), 
                  some(colnames(object@perSampleMetadata[[i]])), "\n") 
            }
          })

setValidity("MultiAssayExperiment", 
            function(object) {
              msg <- validMsg(NULL, NULL)
              perSampleOk <- lapply(lapply(object@perSampleMetadata, colnames),
                                    identical, y=colnames(object))
              if (any(!perSampleOk)) {
                offenders <- names(perSampleOk)[!perSampleOk]
                msg <- validMsg(msg,
                                paste("Invalid perSampleMetadata elements:",
                                      paste(offenders, collapse=", ")))
                msg <- validMsg(msg, 
                                paste("colnames(x@perSampleMetadata)$item must",
                                      "match colnames(x) for each & every item",
                                      "present in x@perSampleMetadata."))
              }
              if (is.null(msg)) TRUE else msg
            })

##  This doesn't work yet.  In theory, it should, but we all know how that goes.
##
## setMethod("[", c("MultiAssayExperiment", "ANY", "ANY"),
##        function(x, i, j, ..., drop=TRUE) {
##          y <- as(x, "RangedSummarizedExperiment")[i, j, ..., drop=drop]
##          y@perSampleMetadata <- lapply(x@perSampleMetadata, function(z) z[,j])
##          as(y, "MultiAssayExperiment")
##        })
##

setAs("RangedSummarizedExperiment", "MultiAssayExperiment", 
      function(from) {
        to <- GenomicRanges:::clone(from)
        class(to) <- "MultiAssayExperiment"
        return(to)
      })

convertToMultiAssay <- function(x) { 
  stopifnot(is(x, "RangedSummarizedExperiment"))
  x <- as(x, "MultiAssayExperiment")
  x@perSampleMetadata <- lapply(x@metadata, function(y) y[, colnames(x)])
  x@metadata <- list()
  return(x)
}

## from AOCS example: this is ugly and should be finessed 
# AOCS_multi <- convertToMultiAssay(AOCS)
