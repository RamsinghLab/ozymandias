#' A master container for linked experiments.  The idea is that the parent has 
#' no assays, and the children have no covariates.  So initializing one of these
#' would usually occur by "lifting" the colData out of a SummarizedExperiment
#' and into the LinkedExperiments, then "demoting" everything else about it 
#' into the child containers.  I think the validity method will end up checking
#' that assays() is empty (heh) and that all child assays have colnames like 
#' the parent's.  This could be relaxed somehow, but I'm not sure how just yet.
#' 
#' @slot linkedAssays   multiple assays, presumably for each sample (or add NAs)
#' 
#' @import    SummarizedExperiment
#' 
#' @export
#'
setClass("LinkedExperiment", 
         contains="SummarizedExperiment0",
         representation(linkedAssays="list"))

setValidity("LinkedExperiment", 
            function(object) {
              msg <- validMsg(NULL, NULL)

              assaysEmpty <- length(assays(object)) == 0 
              if (!assaysEmpty) 
                msg <- validMsg("The assays slot should be empty.")

              colnamesOk <- lapply(lapply(object@linkedAssays, colnames),
                                   identical, y=colnames(object))
              if (any(!colnamesOk)) {
                offenders <- names(colnamesOk)[!colnamesOk]
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
