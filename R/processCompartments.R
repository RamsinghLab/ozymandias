## automate compartment prediction by Sample_Group 
processCompartments <- function(x, exclude=c("chrX", "chrY"), ...) {
  chroms <- seqlevels(x)
  names(chroms) <- chroms
  chroms <- setdiff(chroms, exclude)
  comps <- GRangesList(lapply(chroms, function(y) compartments(x, chr=y)))
  metadata(x)$compartments <- comps
  invisible(x)
}
