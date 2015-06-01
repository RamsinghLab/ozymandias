getPWMenrichment <- function(symbolOrGR, 
                             motifName=NULL,
                             upstream=1500,
                             downstream=200,
                             model=c("gene","transcript"),
                             assembly=c("hg19","mm9")) { 

  model <- match.arg(model)
  if (model != "gene") stop("Transcript mode is not yet supported.")
  assembly <- match.arg(assembly)
  if (assembly == "hg19") { 
    # {{{ human
    library(BSgenome.Hsapiens.UCSC.hg19)
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    library(PWMEnrich.Hsapiens.background)
    data(PWMLogn.hg19.MotifDb.Hsap)
    motifs <- PWMLogn.hg19.MotifDb.Hsap
    library(org.Hs.eg.db)
    egids <- org.Hs.egSYMBOL2EG
    BSg <- Hsapiens 
    # }}}
  } else if (assembly == "mm9") { 
    stop("Only hg19 is supported at the moment.")
    ## {{{ mouse, eventually, but not yet 
    library(BSgenome.Mmusculus.UCSC.mm9)
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    library(PWMEnrich.Mmusculus.background)
    data(PWMLogn.mm9.MotifDb.Mmus)
    motifs <- PWMLogn.mm9.MotifDb.Mmus
    library(org.Mm.eg.db)
    egids <- org.Mm.egSYMBOL2EG
    BSg <- Mmusculus
    # }}}
  } else { 
    stop(paste(assembly, "is not supported at the present time"))
  }

  ## get the coordinates of the promoter and its reference sequence 
  if (!is(symbolOrGR, "GenomicRanges")) { 
    symbol <- symbolOrGR
    egid <- get(symbol, egids)
    overall <- upstream + downstream
    txs <- transcriptsBy(txdb, model)[[egid]]
    promoter <- reduce(resize(flank(txs, upstream), overall, fix="start"))
  } else { 
    ##
    ## e.g.
    ##
    ## c22orf26dmr <- GRanges('chr22', IRanges(start=46400000, end=46515000))
    ## getPWMenrichment(c22orf26dmr)
    ## 
    stopifnot(length(symbolOrGR) == 1) 
    promoter <- symbolOrGR
    symbol <- paste0(seqnames(symbolOrGR), ":", 
                     start(symbolOrGR), "-", 
                     end(symbolOrGR))
  }
  baseseq <- getSeq(BSg, promoter)

  if(is.null(motifName)) {
    ## see what motifs are most highly enriched in the promoter sequence 
    res <- motifEnrichment(baseseq, PWMLogn.hg19.MotifDb.Hsap)
    promoter_report <- groupReport(res)
    targetedMotifs <- which(!is.na(promoter_report$target))[1:20] 
    plot(promoter_report[targetedMotifs], fontsize=7, id.fontsize=10)
    invisible(promoter_report)
  } else { 
    ## look for binding sites of a specific factor/motif in the promoter 
    id <- grep(motifName, ignore=T, names(motifs$pwms), value=T)
    pwms <- motifs$pwms[id]
    scores <- motifScores(baseseq, pwms, raw.scores=T)
    plotMotifScores(scores, cols=c("green", "red"), 
                    main=paste(symbol, "promoter motif scores for", motifName))
    invisible(scores)
  }
}
