#' grab the mitochondrial reads from a BAM and estimate their fraction
#' nb. this could probably be done faster for a list of BAMs but it's not
#'
#' @param bam     a BAM filename, which should have an index (else index it)
#' @param chrM    what the mitochondrial contig is called. Default is "chrM" 
#'
#' @import GenomicAlignments
#' @import Rsamtools
#'
#' @export
getMT <- function(bam, chrM="chrM") {

  bai <- paste0(bam, ".bai")
  if (!file.exists(bai)) indexBam(bam)
  bamfile <- BamFile(bam, index=bai, asMates=TRUE)
  idxStats <- idxstatsBam(bamfile)
  rownames(idxStats) <- idxStats$seqnames
  mtFrac <- idxStats[chrM, "mapped"] / sum(idxStats[, "mapped"])
  mtRange <- GRanges(chrM, IRanges(1, idxStats[chrM, "seqlength"]), "*")
  mtView <- BamViews(bam, bai, bamRanges=mtRange)
  mtParam <- ScanBamParam(what=c("seq","mate_status"))
  mtReads <- readGAlignments(mtView, param=mtParam)[[1]]
  attr(mtReads, "mtFrac") <- mtFrac
  return(mtReads)

}
