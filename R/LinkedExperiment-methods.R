#' @describeIn LinkedExperiment
#' 
#' A referential-integrity-preserving subset operator
#'
#' @param x       a LinkedExperiment
#' @param i       optional row index or indices (will be ignored)
#' @param j       optional column index or indices (will be used) 
#'
#' @return        A valid LinkedExperiment w/column-matching linked assays
#' 
#' @details       If and only if 'j' is set, the LinkedExperiment 'x' 
#'                briefly becomes a 'RangedSummarizedExperiment' while its 
#'                linkedAssays are recursively subsetted. The two pieces are 
#'                then put back together as a valid LinkedExperiment.
#' 
#' @export
#'
setMethod("[", c("LinkedExperiment", "ANY", "ANY"),
          function(x, i, j, ..., drop=TRUE) {
            ## Only if linkedAssays isn't empty and j is set
            if (!missing(j) && length(x@linkedAssays) > 0) {
              x@linkedAssays <- lapply(x@linkedAssays, 
                                       function(z) z[,j])
            }
            callNextMethod()
          })

# FIXME: steal SummarizedExperiment's "show" and drop the assays line
setMethod("show", "LinkedExperiment", 
          function(object) { # {{{
            callNextMethod()
            some <- function(x) { 
              ifelse(length(x) %in% 1:3, 
                     paste(x, collapse=", "),
                     paste(x[1], "...", x[length(x)], sep=", "))
            }
            cat("linkedAssays:", names(object@linkedAssays), "\n")
            for (i in names(object@linkedAssays)) { 
              cat(paste0("  colnames(", i, "):"), 
                  some(colnames(object@linkedAssays[[i]])), "\n") 
            } 
          } # }}}
)

## will eventually rename linkedAssays slot to match this method 
setGeneric("linkedAssays", function(x) standardGeneric("linkedAssays"))
setGeneric("linkedAssays<-", function(x, y) standardGeneric("linkedAssays<-"))
setMethod("linkedAssays", "LinkedExperiment", 
          function(x) x@linkedAssays)

setReplaceMethod("linkedAssays", c("LinkedExperiment", "ANY"), 
          function(x, y) { 
            if (identical(colnames(x), colnames(y))) { 
              x@linkedAssays <- append(x@linkedAssays, y)
            } else { 
              stop("Column names do not match!")
            }
            x
          })

setAs("RangedSummarizedExperiment", "LinkedExperiment", 
      function(from) {
        new("LinkedExperiment", 
            SummarizedExperiment(list(), 
                                 colData=colData(from),
                                 metadata=metadata(from)),
            list(RENAME_ME=SummarizedExperiment(assays(from),
                                                rowRanges=rowRanges(from))))
      })
