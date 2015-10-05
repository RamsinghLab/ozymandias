#' as the name says, if there are IDATs attached to a grSet, rerun them.
#'
#' @param grset   a grSet with $supplementary_file 
#' @param name    the name of the resulting rgSet and grSet files
#'
#' @return        a reanalyzed and funnorm()'ed grSet
#' 
#' @export
#'
reanalyzeFromIDATs <- function(grSet, name) {
  grSet <- prepForReanalysis(grSet)
  if (!is.null(grSet$Basename)) {
    targets <- as(colData(grSet), "data.frame")
    rgSet <- processIDATs(targets, name=name)
    newGrSet <- processMeth(rgSet, name=name)
  }
  return(newGrSet)
}
