#' grab the mitochondrial reads from a BAM and estimate their fraction
#' nb. this could probably be done faster for a list of BAMs but it's not
#' nb. nb. returns NuMt-depleted mitochondrial GenomicAlignments
#'
#' @param bam       a BAM filename, which should have an index (else index it)
#' @param chrM      what the mitochondrial contig is called. Default is "chrM" 
#' @param mtGenome  what mitochondrial assembly was used (default is hg19) 
#' @param plotMAPQ  plot distribution of mitochondrial mapping quality? (FALSE)
#'
#' @import GenomicAlignments
#' @import Rsamtools
#'
#' @export
getMT <- function(bam, chrM="chrM", mtGenome="hg19", plotMAPQ=FALSE) {

  bai <- paste0(bam, ".bai")
  if (!file.exists(bai)) indexBam(bam)
  bamfile <- BamFile(bam, index=bai, asMates=TRUE)
  idxStats <- idxstatsBam(bamfile)
  rownames(idxStats) <- idxStats$seqnames
  mtFrac <- idxStats[chrM, "mapped"] / sum(idxStats[, "mapped"])
  message(bam, " has ~", round(mtFrac * 100, 2), "% mitochondrial reads.")

  mtRange <- GRanges(chrM, IRanges(1, idxStats[chrM, "seqlength"]), "*")
  mtView <- BamViews(bam, bai, bamRanges=mtRange)
  if (!base::grepl(mtGenome, bam)) {
    message(mtGenome, " (supplied or default) isn't found in your bam filename")
  }
  flags <- scanBamFlag(isPaired=TRUE, 
                       isProperPair=TRUE, 
                       isUnmappedQuery=FALSE, 
                       hasUnmappedMate=FALSE, 
                       isSecondaryAlignment=FALSE, 
                       isNotPassingQualityControls=FALSE, 
                       isDuplicate=FALSE)

  mtParam <- ScanBamParam(flag=flags, what=c("seq","mapq")) 
  mtReads <- suppressWarnings(readGAlignments(mtView, param=mtParam)[[1]])
  mtReadLength <- median(width(mcols(mtReads)$seq))
  attr(mtReads, "mtReadLength") <- mtReadLength
  attr(mtReads, "mtFrac") <- mtFrac
  genome(mtReads) <- mtGenome
  mtReads <- keepSeqlevels(mtReads, chrM)
  isCircular(seqinfo(mtReads)) <- TRUE 
  seqinfo(bamRanges(mtView)) <- seqinfo(mtReads)[seqlevels(bamRanges(mtView))]
  attr(mtReads, "mtView") <- mtView

  if (plotMAPQ) {
    plot(density(mcols(mtReads)$mapq), type="h", col="red",
         xlab="MAPQ", ylab="fraction of reads with this MAPQ", 
         main=paste("Mitochondrial read mapping quality for\n", bam))
  }

  return(mtReads)

}
