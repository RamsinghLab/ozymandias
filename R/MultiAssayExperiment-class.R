#' Something like a SummarizedExperiment, but with referential integrity checks.
#' This is mostly a proof of concept for a well-behaved integration container.
#' 
#' @slot perSampleMetadata associated per-sample data, checked against colnames
#' 
#' @import    SummarizedExperiment
#' 
#' @details   perSampleMetadata will shortly be renamed to "linkedAssays".
#'
#' @export
#'
setClass("MultiAssayExperiment", 
         contains="RangedSummarizedExperiment", # finesse the RSE issue?
         representation(perSampleMetadata="list")) # rename to linkedAssays 

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
            }
)
