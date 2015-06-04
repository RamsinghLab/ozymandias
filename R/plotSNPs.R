#' plot SNP probes from DNA methylation arrays to avoid label swaps
#' 
#' @param x           SummarizedExperiment-like object with DNA methylation data
#' @param individuals How many distinct individuals are (allegedly) in the data?
#' @param rotate      rotate the plot? (default is FALSE) 
#' 
#' @return            subject identity assignments for the dataset 
#'
plotSNPs <- function(x, individuals=NULL, rotate=FALSE, ...) {

  name <- as.character(match.call()["x"])

  tmp <- matrix()
  if (is(x, "RangedSummarizedExperiment")) { # {{{ sanity checks
    probes <- grep("^rs", names(SummarizedExperiment::rowRanges(x)))
    snps <- grep("^rs", names(SummarizedExperiment::rowRanges(x)), val=T)
    tmp <- assays(x, withDimnames=F)[[1]][probes, ]
    if (nrow(tmp) < 2 && "SNPs" %in% names(metadata(x))) {
      if (!identical(colnames(metadata(x)$SNPs), colnames(x))) {
        if (all(colnames(x) %in% colnames(metadata(x)$SNPs))) {
          message("Additional columns in SNP matrix not found in SE, dropping.")
        } else { 
          stop("Colnames of your SNP matrix do not match those of your data!")
        }
      }
      tmp <- metadata(x)$SNPs[, colnames(x)] 
    }
  } # }}}
  if (nrow(tmp) < 2) {
    stop("Need a SummarizedExperiment with SNP (rsXX) features... none found")
  }

  ## fit w/mclust
  SNPfit <- Mclust(as.numeric(tmp), G=3)
  calls <- matrix(SNPfit$class - 1, ncol=ncol(tmp))
  rownames(calls) <- rownames(tmp) 
  colnames(calls) <- colnames(tmp) 

  if (is.null(individuals)) individuals <- ncol(x) / 2
  SNP <- c("blue","yellow","red")
  heading <- paste("SNPs for", individuals, "individuals in", name)
  message("Tracking plot for clusters...")
  rcc <- ConsensusClusterPlus(calls, maxK=individuals, reps=100, tmyPal=SNP,
                              distance="manhattan", clusterAlg="pam", 
                              ...)[[individuals]]
  individual <- as.factor(rcc$consensusClass)
  par(mfrow=c(1,1)) 
  if(rotate) calls <- t(calls)
  heatmap(calls, col=SNP, scale="none", main=heading)
  message("Assigned identities:")
  print(table(individual))
  invisible(individual)
}
