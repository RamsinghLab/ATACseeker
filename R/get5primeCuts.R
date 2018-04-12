#' get the 5' cut sites from a GAlignmentPairs (usually of an ATAC BAM) 
#'
#' @param   galp        a GAlignmentPairs object
#' @param   maxSize     maximum insert size (for fragmentLengths(galp))
#' @param   shrink      move the reads around like in the original paper? (TRUE)
#' @param   shrinkBy    shrink by how much? (realistically this is pointless...)
#'
#' @return  the 5' ends, from getEnds() 
#'
#' @examples
#' library(Homo.sapiens)
#' CUX1_EGID <- org.Hs.egSYMBOL2EG[['CUX1']]
#' CUX1 <- reduce(transcriptsBy(Homo.sapiens, 'gene')[[CUX1_EGID]])
#' CUX1bump <- GRanges('chr7', IRanges(101499132, 101501052), '+')
#'
#' @import GenomicAlignments
#' 
#' @export
get5primeCuts <- function(galp, maxSize=1000, shrink=TRUE, shrinkBy=NULL, ...) {

  if(!is.null(maxSize)) {
    whichPairs <- which(fragmentLengths(galp, plotAs='none') <= maxSize)
    galp <- galp[whichPairs] 
  }

  ## median of the fragment length distribution
  if(is.null(shrinkBy)) shrinkBy <- round(median(width(granges(galp))))

  ## get the actual cuts
  if(shrink == TRUE) {
    res <- getEnds(shrinkPairs(granges(galp), shrinkBy))
  } else {
    res <- getEnds(granges(galp))
  }
  attributes(res)$metadata$what <- 
    paste("5' cuts,", ifelse(shrink,"shrunken","unshrunken"))
  return(res)

}
