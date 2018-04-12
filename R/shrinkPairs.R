#' This seems to work better than anything else for "shifting" cuts
#'
#' @param x         a GAlignmentPairs or GRanges
#' @param shrinkBy  shrink by this many bases
#' 
#' @return          the input object, widths shrunken by shrinkBy bases
#'
#' @import GenomicAlignments
#'
#' @export
shrinkPairs <- function(x, shrinkBy) {
  ## could probably do this better with a normexp/fixseq approximation
  resize(x, pmax(width(x) - shrinkBy, shrinkBy), fix='center')
}
