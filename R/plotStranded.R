#' convenience function
#'
#' @param x     coverage
#' @param gr    a single GRange
#' @param ...   args to pass along to plotCoverage
#'
#' @import GenomicRanges
#'
#' @export
plotStranded <- function(x, gr, ...) plotCoverage(x, gr, ..., stranded=TRUE)

