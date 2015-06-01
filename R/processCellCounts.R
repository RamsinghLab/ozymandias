## initially, this is for blood, ONLY.  It could handle brains
## 
processCellCounts <- function(grSet, name, Tissue=NULL, ...) {
  grSet <- processDNAmAge(grSet, name=name, Tissue=Tissue)
  colData(grSet) <- cbind(colData(grSet), DataFrame(getBloodCellCounts(grSet)))
  return(grSet)
}
