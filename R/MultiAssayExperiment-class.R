#' Something like a SummarizedExperiment, but with "holes" and metadata checks
#'
#' @slot perSampleMetadata associated per-sample data, checked against colnames
#'
#' @export
#'
setClass("MultiAssayExperiment", 
         representation(perSampleMetadata="List"), # or a DataFrameList maybe?
         contains="RangedSummarizedExperiment") # RangedSummarizedExperiments?
setMethod("show", "MultiAssayExperiment", 
          function(object) {
            callNextMethod()
            cat("perSampleMetadata:", names(object@perSampleMetadata), "\n")
          })
setValidity("MultiAssayExperiment", 
            function(object) {
              msg <- validMsg(NULL, NULL)
              perSampleOk <- lapply(lapply(object@perSampleMetadata, colnames),
                                    identical, y=colnames(object))
              if (any(!perSampleOk)) {
                offenders <- names(perSample)[!perSample]
                msg <- validMsg(msg,
                                paste("Invalid perSampleMetadata elements:",
                                      paste(offenders, collapse=", ")))
              }
              if (is.null(msg)) TRUE else msg
            })

## FIXME: add subsetting methods that ensure validity 

# setAs("RangedSummarizedExperiment", "MultiAssayExperiment", 
#      function(from) callNextMethod()) ??

## from AOCS example: this is ugly and should be finessed 
# AOCS_multi <- as(AOCS, "MultiAssayExperiment")
# AOCS_multi@perSampleMetadata <- List(AOCS@metadata)
# AOCS_multi@metadata$miRNA <- NULL
