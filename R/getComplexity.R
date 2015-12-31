#' complexity estimates from a GAlignmentPairs object, formatted for preseqR 
#'
#' @param x   the GAlp object or GRanges
#' 
#' @return  a matrix of estimates
#' 
#' @export
getComplexity <- function(x) {

  if (is(x, 'GAlignmentPairs')) x <- getFeatureCounts(as(x, 'GRanges'))
  if (!'seen' %in% names(x)) x <- getFeatureCounts(x)
  freq <- seq_len(max(range(x$seen)))
  breaks <- c(0, freq)
  d <- data.frame(freq=freq, species=hist(x$seen, breaks=breaks, plot=F)$counts)
#  d <- d[ which(d[,2] != 0), ] 
  return(data.matrix(d))

}
