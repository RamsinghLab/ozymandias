#' @describeIn MultiAssayExperiment
#' 
#' A referential-integrity-preserving subset operator
#'
#' @param x       a MultiAssayExperiment
#' @param i       optional index or indices
#' @param j       optional index or indices
#'
#' @return        A valid MultiAssayExperiment w/column-matching linked assays
#' 
#' @details       If and only if 'j' is set, the MultiAssayExperiment 'x' 
#'                briefly becomes a 'RangedSummarizedExperiment' while its 
#'                linkedAssays are recursively subsetted. The two pieces are 
#'                then put back together as a valid MultiAssayExperiment.
#' 
#' @export
#'
setMethod("[", c("MultiAssayExperiment", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE) {
            ## Only if perSampleMetadata isn't empty and j is set
            if (!missing(j) && length(x@perSampleMetadata) > 0) {
              x@perSampleMetadata <- lapply(x@perSampleMetadata, 
                                            function(z) z[,j])
            }
            callNextMethod()
          })

setMethod("show", "MultiAssayExperiment", 
          function(object) { # {{{
            callNextMethod()
            some <- function(x) { 
              ifelse(length(x) %in% 1:3, 
                     paste(x, collapse=", "),
                     paste(x[1], "...", x[length(x)], sep=", "))
            }
            cat("linkedAssays:", names(object@perSampleMetadata), "\n")
            for (i in names(object@perSampleMetadata)) { 
              cat(paste0("  colnames(", i, "):"), 
                  some(colnames(object@perSampleMetadata[[i]])), "\n") 
            } 
          } # }}}
)

## will eventually rename perSampleMetadata slot to match this method 
setGeneric("linkedAssays", function(x) standardGeneric("linkedAssays"))
setGeneric("linkedAssays<-", function(x, y) standardGeneric("linkedAssays<-"))

setMethod("linkedAssays", "MultiAssayExperiment", 
          function(x) x@perSampleMetadata)

setReplaceMethod("linkedAssays", c("MultiAssayExperiment", "ANY"), 
          function(x, y) { 
            if (identical(colnames(x), colnames(y))) { 
              x@perSampleMetadata <- append(x@perSampleMetadata, y)
            } else { 
              stop("Column names do not match!")
            }
            x
          })

## phase these out 
setGeneric("perSampleMetadata", 
           function(object) standardGeneric("perSampleMetadata"))
setMethod("perSampleMetadata", "MultiAssayExperiment", 
           function(object) object@perSampleMetadata)

setAs("RangedSummarizedExperiment", "MultiAssayExperiment", 
      function(from) new("MultiAssayExperiment", from, list()))
