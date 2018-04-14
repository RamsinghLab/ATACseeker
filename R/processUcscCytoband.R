#' utility function to pull in UCSC cytoband / arm files 
#' this is then applied to the mcols() of a GRanges or GAlignmentPairs
#' and then one can trivially parallelize over cytobands or arms with split() 
#' 
#' @param   filename    the UCSC cytoband file
#'
#' @return  a GRanges with $chrArm and $chrBand as mcols
#' 
#' @import GenomicRanges 
#' 
#' @export
processUcscCytoband <- function(filename) {
  cb <- read.table(filename, 
                   col.names=c('chrom','chromStart','chromEnd','band','stain'))
  cb$chrArm <- paste0(cb$chrom, substr(cb$band, 1, 1))
  cb$chrBand <- paste0(cb$chrom, cb$band)
  makeGRangesFromDataFrame(cb, keep.extra.columns=TRUE)
}
