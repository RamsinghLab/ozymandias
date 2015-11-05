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

  ## the following is obviated by MultiAssayExperiment
  if (class(x) != "MultiAssayExperiment") {
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
  } else {
    tmp <- perSampleMetadata(x)$SNPs 
  }

  if (nrow(tmp) < 2) {
    stop("Need a SummarizedExperiment with SNP (rsXX) features... none found")
  }


  # fit w/mclust
  SNPfit <- Mclust(as.numeric(tmp), G=3)
  calls <- matrix(SNPfit$class - 1, ncol=ncol(tmp))
  rownames(calls) <- rownames(tmp) 
  colnames(calls) <- colnames(tmp) 

  SNP <- c("blue","yellow","red")
  heading <- paste("SNPs for", ncol(x) , "samples in", name)

  par(mfrow=c(1,1)) 
  if(rotate) calls <- t(calls)
  if (!is.null(individuals) || (ncol(x) < 99)) { 
    if (is.null(individuals)) individuals <- ncol(x) / 2
    message("Tracking plot for clusters...")
    rcc <- ConsensusClusterPlus(calls, maxK=individuals, reps=10, tmyPal=SNP,
                                distance="manhattan", clusterAlg="pam", 
                                ...)[[individuals]]
    individual <- as.factor(rcc$consensusClass)
    individuals <- length(levels(individual))
    heading <- paste("SNPs for", ncol(x) , "samples from", 
                     individuals, "subjects in", name)
  }

  heatmap(calls, col=SNP, scale="none", main=heading)
  
  if (!is.null(individuals)) {
    message("Assigned identities:")
    print(table(individual))
    invisible(individual)
  } else { 
    invisible(1:ncol(x))
  }

}
