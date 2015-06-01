## prep beta values, annotations for uploading to Steve Horvath's calculator 
## eventually I will wrap the damned thing properly to do it internally
## 
processDNAmAge <- function(grSet, name, Tissue=NULL, ...) {
 
  cols <- c("Sample_Name", "Age", "Female", "Tissue")

  if ("predictedSex" %in% names(colData(grSet)) &&  # {{{
      !"Female" %in% names(colData(grSet))) {
    grSet$Female <- as.numeric(grepl("^F", ignore=TRUE, grSet$predictedSex))
  } # }}}
  if (!is.null(Tissue) && !"Tissue" %in% names(colData(grSet))) { # {{{
    grSet$Tissue <- Tissue
  } # }}}
  if (!"Age" %in% names(colData(grSet))) { # {{{
    grSet$Age <- NA
  } # }}}
  if (any(grepl("HumanMethylation450", annotation(grSet)))) { # {{{
    data(datMini450k) ## from Horvath, subset of datMiniAnnotation.csv
    loci <- rownames(datMini450k) # }}}
  } else if (any(grepl("HumanMethylation27", annotation(grSet)))) { # {{{
    data(datMini27k) ## from Horvath, subset of datMiniAnnotation.csv
    loci <- rownames(datMini27k) # }}}
  } else { # {{{ cannot proceed, just tidy up and return the grSet 
    message("Cannot dump appropriate data for DNAm age predictor.")
    invisible(grSet) ## return a tidied version of the original
    message("Returning a (partially-reannotated) GenomicRatioSet.")
  } # }}}

  if (all(cols %in% names(colData(grSet))) && all(loci %in% rownames(grSet))) {
    annotFilename <- paste0(name, "_DNAmAge_Annotations.csv")
    dnamageFilename <- paste0(name, "_DNAmAge_BetaValues.csv") 
    message("Writing ", length(loci), " loci to ", dnamageFilename, "...")
    write.csv(getBeta(grSet[loci, ]), dnamageFilename)
    message("Writing annotations to ", annotFilename, "...")
    write.csv(colData(grSet)[, cols], annotFilename)
    message("Done.")
    invisible(grSet) ## return a tidied version of the original
  } else {
    message("Cannot dump appropriate data for DNAm age predictor.")
    message("Returning a (partially-reannotated) GenomicRatioSet.")
    invisible(grSet) ## return a tidied version of the original
  }

}
