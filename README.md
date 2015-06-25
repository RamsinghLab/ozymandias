# ozymandias

There are a lot of obvious QC/label matching steps missing from currently available/visible DNA methylation microarray pipelines, which while not being as fashionable or expensive as single-cell BSseq etc., benefit from tens of thousands of existing samples in hundreds of existing experiments.  After a few dozen kicks and stings I ended up writing a package to handle batches of DNA methylation measurements (which, if ascertained from Illumina arrays, are accompanied by copy number information and high-MAF SNPs), possibly accompanied by matched RNAseq or mRNA microarray measurements, translocation/inversion/mutation covariates, and tissue type indicators.  

The goal here is a preprocessing step that makes it much harder to screw things up.  So, for example, if the chrX/chrY copy number doesn't agree with the X inactivation status by DNA methylation, that gets flagged.  If the high-MAF SNPs on an array don't agree across samples supposedly from the same person, that gets flagged.  If the "epigenetic age" of normal samples is wildly different from their specified age, that gets flagged.  If copy number aberrations for a supposedly identical sample are radically different between the RNA and DNA samples, that gets flagged.  This is all terribly boring, which might explain why many groups seem not to bother with it.  That, in turn, might have something to do with the difficulty reported in reproducing flashy findings in glamor journals...  

If you've never been burned by label swaps or horrible batch effects or human error, go right ahead and ignore all of these things.  What could possibly go wrong?  Other than losing a year or two of your life chasing ghosts, that is...
