#' function to fit Mas-o-Menos (MoM) predictor scores (usually) against DMRs
#' 
#' @param atDMRs    SummarizedExperiment-derived object collapsed across regions
#' @param tumorName Optional tumor name for the plots (default: guess from args)
#' 
#' @return          a list with K-M and Cox fits, effect sizes & mm coefficients
#'
fitMoM <- function(atDMRs, tumorName=NULL) {

  if(is.null(tumorName)) tumorName <- as.character(match.call()["atDMRs"])

  atDMRs <- atDMRs[, hasOS(atDMRs)] 
  effectSizes <- rowCoxTests(assays(atDMRs)$Beta, OSurv(atDMRs), option="slow")
  mm.coefs <- masomenos(t(assays(atDMRs)$Beta), OSurv(atDMRs), option="slow")
  atDMRs$mm.score <- t(assays(atDMRs)$Beta) %*% mm.coefs
  mm.fit <- Mclust(atDMRs$mm.score)
  plot(density(atDMRs$mm.score), xlab="score", ylab="frequency", 
       main=paste("risk score distribution for", name))
  readline("Press enter to continue...")
  if(length(unique(mm.fit$class)) == 1) {
    message("Splitting on median...")
    med <- median(atDMRs$mm.score)
    atDMRs$mm.class <- ifelse(atDMRs$mm.score > med, 2, 1)
  } else {
    message("Using mclust to assign groups...")
    plot(density(atDMRs$mm.score),
         xlab="methRisk score",
         main=paste("Clustering of DNAm scores:", tumorName))
    atDMRs$mm.class <- mm.fit$class
  } 

  fit <- survfit(OSurv(atDMRs) ~ atDMRs$mm.class)

  ## add option to fit DNAmAge as well?
  if ("age" %in% names(colData(atDMRs))) {
    age <- TRUE
    cph <- coxph(OSurv(atDMRs) ~ atDMRs$mm.class + atDMRs$age)
  } else { 
    age <- FALSE 
    cph <- coxph(OSurv(atDMRs) ~ atDMRs$mm.class)
  }
  p <- round(summary(cph)$coefficients["atDMRs$mm.class", 5], 3)
  pTitle <- ifelse(age, paste("age-adjusted p =", p), paste("p =", p))
  plot(fit, col=c("green", "red"), lwd=3,
            main=paste(tumorName, "survival, by methRisk category:", pTitle),
            ylab="Fraction surviving",
            xlab="Days surviving")
  return(list(fit=fit, cph=cph, mm.coefs=mm.coefs, effectSizes=effectSizes))

}
