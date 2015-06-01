## I think Vince and Levi and pals were working on a MergedExperiment thingy...
## It would be a good idea to find out how far that has got and/or use it here.
##
## mandatory columns: target$Sample_Group; target$Sample_Name; target$Basename 
## what they contain: "control" vs. others; title of sample; barcode of sample
## optional columns: whatever other covariates each sample may have available!
## 
get450kAndCNV <- function(targets, DMRs=TRUE, blood=FALSE, CNV=FALSE, 
                          RNA=c(FALSE, "none", "affy", "rnaseq"), 
                          ctl="control", ...) {

  name <- as.character(match.call()["targets"])
  grSetName <- paste0(name, "_450k.rds")
  message("Please note that the name of the variable holding your targets,")
  message(name)
  message("will we used as the base name for all output files, e.g.")
  message(grSetName)
  message("If this is a problem, ABORT THE FUNCTION NOW and rename.")
  ## no? alright then...  (maybe add a readline("Press enter to proceed:") here)
  message("\n\nProcessing ", name, "...")

  RNA <- match.arg(RNA)
  cols <- c("Sample_Group","Sample_Name","Basename")
  if (RNA %in% c("affy","rnaseq")) cols <- c(cols, c("RNAfile","RNAplatform"))
  if (!all(cols %in% names(targets))) { 
    # {{{ halt and catch fire 
    message("Your targets data.frame is missing mandatory columns:")
    for (missed in setdiff(cols, names(targets))) message(missed)
    stop("Please fix this and then try again...")
    # }}}
  } else { 
    for (acol in cols) targets[, acol] <- as.character(targets[, acol])
  }

  rgSet <- processIDATs(targets, name) ## pulls in SNPs and detection p-values
  grSet <- processMeth(rgSet, name) ## masks betas by pval > 0.01 & keeps SNPs

  ## the most time-consuming step:
  if (CNV == TRUE) { 
    metadata(grSet)$CNVs <- processCNV(rgSet, name)
    saveRDS(grSet, grSetName)
  }

  ## add RNAseq or Affy mRNA/lncRNA?
  if (RNA %in% c("affy", "rnaseq")) { # {{{
    if (RNA == "affy") { 
      ## currently fRMA, would like to switch to SCAN.UPC
      metadata(grSet)$affy <- processAffy(grSet$RNAfile)
      saveRDS(grSet, grSetName)
    } else if (RNA == "rnaseq") { 
      if (length(unique(grSet$RNAfile) > 1) || 
          length(unique(grSet$RNAplatform) > 1)) {
        stop("You have > 1 RNAseq platforms listed; this will not work.")
      } else {
        RNAfile <- unique(grSet$RNAfile) ## transcript counts
        RNAplatform <- unique(grSet$RNAplatform) ## transcript annotations
        metadata(grSet)$rnaseq <- processRNAseq(RNAfile, RNAplatform)
      }
    }
  } # }}}

  ## optional D/VMRs 
  if (DMRs == TRUE) { 
    ## {{{ call D/VMRs 
    design <- NULL
    ctl <- "control"
    if (any(grSet$Sample_Group == ctl) && !all(grSet$Sample_Group == ctl)) {
      cols <- c("dmrCase", "Female")  ## {{{ call simple DMRs 
      grSet$dmrCase <- as.numeric(grepl("control", grSet$Sample_Group)) # }}}
      if (!"Female" %in% names(colData(grSet))) { #{{{
        grSet$Female <- as.numeric(grepl("^F", ignore=T, grSet$predictedSex))
      } # }}}
      if (blood == TRUE) { # {{{ process in a more sophisticated manner
        message("Adjusting for blood cell composition...")
        grSet <- processCellCounts(grSet, name=name, Tissue="Blood PBMC")
        cols <- c(cols, "Age", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran") 
        if (!all(cols %in% names(colData(grSet)))) {
          message("Missing columns in colData, falling back to default design")
        } else { 
          design <- with(as(colData(grSet)[,cols],"data.frame"), 
                         model.matrix(~.))
        }
        message("Done. You may also estimate DNAm age on Horvath's website.")
        message("Try adding")
        message(" ~ DNAmAge + Age + CD8.naive + CD8pCD28nCD45RAn + Female +")
        message("   PlasmaBlast + CD4T + NK + Mono + Gran")
        message("to your design matrix if you want to adjust for these factors")
        # }}}
      } else { # {{{ just call cases vs. controls 
        design <- with(as(colData(grSet)[,cols],"data.frame"), model.matrix(~.))
      } # }}}
    }
    # }}}
    metadata(grSet)$DMRs <- processDMRs(rgSet, design, pcutoff=.1, 
                                        betacutoff=.1, parallel=T)
    metadata(grSet)$VMRs <- processVMRs(rgSet, pcutoff=.1, parallel=T)
    saveRDS(grSet, grSetName)
  } 

  invisible(grSet) 

}
