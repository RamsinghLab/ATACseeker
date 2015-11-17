#' collect various cross-correlation and window size stats for QC
#' 
#' @param curbam    the name of the BAM file to do this for 
#' @param pe.param  paired end parameters, derived elsewhere
#' 
#' @return          the output of csaw::profileSites after massaging inputs
#'
#' @export
collectStats <- function(curbam, pe.param) { 
  wcnt <- windowCounts(curbam, spacing=50, width=50, param=pe.param, filter=20)
  rsums <- rowSums(assays(wcnt)$counts)
  maxed <- findMaxima(rowRanges(wcnt), range=1000, metric=rsums)
  rmax <- rowRanges(wcnt)[maxed]
  weight <- 1 / rsums[maxed]
  profileSites(curbam, rmax, param=pe.param, weight=weight)
}
