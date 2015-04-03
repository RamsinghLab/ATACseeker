## random crap, e.g. processing UCSC cytoband files for easier parallelization
##
processUcscCytoband <- function(filename) {
  cb <- read.table(filename, 
                   col.names=c('chrom','chromStart','chromEnd','band','stain'))
  cb$chrArm <- paste0(cb$chrom, substr(cb$band, 1, 1))
  cb$chrBand <- paste0(cb$chrom, cb$band)
  makeGRangesFromDataFrame(cb, keep.extra.columns=TRUE)
}

byArm <- function(gr) {
  stopifnot('chrArm' %in% names(mcols(gr)))
  reduce(GenomicRanges::split(gr, gr$chrArm))
}

byBand <- function(gr) {
  stopifnot('chrBand' %in% names(mcols(gr)))
  reduce(GenomicRanges::split(gr, gr$chrBand))
}

plusStrand <- function(x) x[which(strand(x) == '+')]
minusStrand <- function(x) x[which(strand(x) == '-')]

grow <- function(x, y) resize(x, width(x) + (y*2))
shrink <- function(x, y) resize(x, pmax(1, width(x) - (y*2)))

