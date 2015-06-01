## ease transition from SummarizedExperiment to RangedSummarizedExperiment
setMethod("update", signature(object = "GenomicRatioSet"),
          function(object) {
              getOldAssay <- function(x, y) slot(x, "assays")[["data"]][[y]]
              grset <- GenomicRatioSet(gr=slot(object, "rowData"),
                                       Beta=getOldAssay(object, "Beta"),
                                       M=getOldAssay(object, "M"),
                                       CN=getOldAssay(object, "CN"),
                                       pData=slot(object, "colData"),
                                       annotation=slot(object, "annotation"),
                                       preprocessMethod=slot(object, "preprocessMethod"))
              metadata(grset) <- as(slot(object, "exptData"), "list")
              return(grset)
          })
