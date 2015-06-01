prepBarcode <- function(x) { 

  xx <- assays(x)$barcode
  if (any(is.na(xx))) {
    message('Imputing NAs...')
    xx <- impute.knn(xx)$data
  }
  rownames(xx) <- rownames(x)
  colnames(xx) <- colnames(x)
  return(xx)
        
}
