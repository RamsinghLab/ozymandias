#' plot X chromosome probes from DNA methylation arrays to confirm gender
#' 
#' @param x       a SummarizedExperiment-like dataset with DNA methylation data
#' @param nFeats  how many features (max) to use for clustering (default 100)
#' @param rotate  rotate the plot? (default is FALSE) 
#' 
#' @return        gender assignments for the subjects 
#'
plotX <- function(x, nFeats=100, rotate=FALSE, ...) {

  if(!(is(x, "SummarizedExperiment") | is(x, "RangedSummarizedExperiment"))) {
    stop("You need to provide a SummarizedExperiment for this to work")
  } else if (!("chrX" %in% runValue(seqnames(x)))) {
    stop("Your SummarizedExperiment does not contain any chrX probes!")
  }

  ## order by SD, then select top nFeats CpG loci 
  chrX <- keepSeqlevels(x, "chrX")
  if (nrow(chrX) < 1) stop("You do not seem to have any ChrX probes!")
  x.chrX.SD <- rowSds(assays(chrX)[[1]]) ## M, beta, who cares
  names(x.chrX.SD) <- rownames(chrX)
  Xloci <- names(head(sort(x.chrX.SD, decreasing=TRUE), nFeats))
  XInd <- which(rownames(x) %in% grep("^cg", Xloci, val=T))

  name <- as.character(match.call()["x"])
  asy <- names(assays(x, withDimnames=F))[[1]]
  tmp <- assays(x, withDimnames=F)[[asy]][XInd, ]

  ## fit w/mclust
  Xfit <- Mclust(t(tmp), G=2)
  gender <- c("M", "F")[Xfit$class]
  meanByGender <- tapply(colMeans(tmp), gender, mean, na.rm=T)
  if (meanByGender["M"] > meanByGender["F"]) {
    ## Flip the gender assignments so they make sense
    gender <- c("F", "M")[Xfit$class]
  }
  sexColors <- c(M="lightblue",F="pink")
  colSide <- sexColors[gender] 
    
  message("Assigned gender:")
  if ("predictedSex" %in% names(colData(x))) {
    print(table(methylationSex=gender, copyNumberSex=colData(x)$predictedSex))
  } else { 
    print(table(gender))
  }

  ## do a test to quantify the chance they"re all the same
  byClust <- split(colMeans(tmp), gender)
  p <- wilcox.test(byClust[["M"]], byClust[["F"]])$p.value
  message("Pr(all samples are the same sex) = ", p)
  jet <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                            "#7FFF7F", "yellow", "#FF7F00", "red", 
                            "#7F0000"))(256)
  par(mfrow=c(1,1)) 
  if(rotate) {
    calls <- t(calls)
    heatmap(tmp, col=jet, scale="none", RowSideColors=colSide,
            main=paste("chrX clustering for", name))
  } else { 
    heatmap(tmp, col=jet, scale="none", ColSideColors=colSide,
            main=paste("chrX clustering for", name))
  }
  invisible(gender)

}
