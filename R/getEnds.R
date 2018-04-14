#' convert a galp into a gr of 5' cuts
#' 
#' @param gr  an intermediate GRanges from a GAlignmentPairs
#'
#' @return    a gr of width-1 cuts
#'
#' @import    GenomicAlignments 
#' 
#' @export
getEnds <- function(gr) { 
  sort(c(resize(gr, 1, fix='start'),
         resize(gr, 1, fix='end')))
}

