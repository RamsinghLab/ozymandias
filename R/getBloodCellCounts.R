## Convenience functions to do "sensible" compositional correction of samples
## These are probably best used in conjunction with Steve Horvath's DNAmAge

## adapted minfi cell count estimation algorithm
## also, fix a bug where it can return counts < 0
## (Houseman's code was less convoluted than Jaffe's in minfi)
getBloodCellCounts <- function(grSet, referenceMset=NULL) { # {{{
  if(is.null(referenceMset)) {
    cat("[estimateCellCounts] Loading reference FlowSorted.Blood.450k...\n")
    data(referenceMset) ## Kere[rownames(FlowSorted.Blood.450k.JaffeModelPars),]
  }
  cellTypes <- c("CD8T","CD4T", "NK","Bcell","Mono","Gran")
  compData <- minfi:::pickCompProbes(referenceMset, cellTypes=cellTypes)
  coefs <- compData$coefEsts
  coefs <- coefs[ intersect(rownames(grSet), rownames(coefs)), ]
  cat("[estimateCellCounts] Estimating composition...\n")
  counts <- minfi:::projectCellType(getBeta(grSet)[rownames(coefs), ], coefs)
  rownames(counts) <- sampleNames(grSet)
  if(any(counts < 0)) {
    counts[ which(counts < 0) ] <- 0
    sums <- rowSums(counts)
    counts <- sweep(counts, 1, sums, "/")
  }
  return(counts) 
} # }}}

## stacked bars (for raw counts)
plotCellCounts <- function(estimates) { # {{{
  require(reshape2)
  if(is(estimates, "matrix")) {
    estimates <- melt(estimates)
    names(estimates) <- c("subject","celltype","fraction")
  } else { 
    stop("Need matrix of counts with columns subject, celltype, fraction")
  }
  require(ggplot2)
  p <- ggplot(estimates, aes(y=fraction, x=factor(subject), fill=celltype)) +
              geom_bar(position="fill", stat="identity") + 
              xlab("subject") + ylab("fraction")
  p <- p + scale_fill_brewer(type="div", palette=7) +
           ggtitle("Estimated leukocyte composition of each sample") +
           theme(panel.background=element_blank(),
                 axis.text.x = element_text(angle=45,hjust=1,vjust=1,size=9))
  p
} # }}}

## e.g.
##
## if(FALSE) {
##   library(minfi)
##   HuGeF.pdat <- readRDS("HuGeF.pData.rds")
##   HuGeF <- read.450k.exp(base=".", targets=HuGeF.pdat)
##   for( i in colnames(HuGeF.counts) ) pData(HuGeF)[,i] <- HuGeF.counts[,i]
##   HuGeF.counts <- estimateCellCounts(HuGeF)
##   HuGeF <- preprocessQuantile(HuGeF)
##   saveRDS(HuGeF, file="HuGeF.rds")
##
##   ## HuGeF age DMRs
##   library(doMC)
##   registerDoMC(2) ## for parallel bump hunting ... be careful though
##   HuGeF.age.DMRs <- bumphunter(HuGeF, dmat, pickCutoff=T, cutoffQ=0.99)
##   saveRDS(HuGeF.age.DMRs, file="HuGeF.age.DMRs.rds")
## }

