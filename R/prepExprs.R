prepExprs <- function(x) { 

  xx <- assays(x)$exprs
  if (any(is.na(xx))) {
    message('Imputing NAs...')
    xx <- impute.knn(xx)$data
  }
  return(xx)
        
}
