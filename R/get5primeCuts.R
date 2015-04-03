##
## library(Homo.sapiens)
## CUX1_EGID <- org.Hs.egSYMBOL2EG[['CUX1']]
## CUX1 <- reduce(transcriptsBy(Homo.sapiens, 'gene')[[CUX1_EGID]])
## CUX1bump <- GRanges('chr7', IRanges(101499132, 101501052), '+')
##
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

##
## This seems to work better than anything else for "shifting" cuts
## 
shrinkPairs <- function(x, shrinkBy) {
  ## could probably do this better with a normexp/fixseq approximation
  resize(x, 
         pmax(width(x) - shrinkBy, shrinkBy), 
         fix='center')
}

##
## convert a galp into a gr of 5' cuts
##
getEnds <- function(gr) { 
  sort(c(resize(gr, 1, fix='start'),
         resize(gr, 1, fix='end')))
}


## exporting 5' cuts after shrinkage:
##
## getDepthAsScore <- function(x) {
##   xx <- x
##   strand(xx) <- '*'
##   xx.uniq <- unique(xx)
##   xx.uniq$score <- countOverlaps(xx.uniq, xx)
##   return(xx.uniq)
## }
## 
## SMO.cuts <- get5primeCuts(galp)
## strand(SMO.cuts) <- '*'
## SMO.cuts.uniq <- unique(SMO.cuts)
## SMO.cuts.uniq$score <- countOverlaps(SMO.cuts.uniq, SMO.cuts)
## export(SMO.cuts.uniq, 'A1.SMO.shrunkenCuts.mm10.wig')
## wigToBigWig('A1.SMO.shrunkenCuts.mm10.wig', seqinfo(Mus.musculus))

get5primeCutDensity <- function(x, smooth=TRUE, sigma=50) {

  ## this seems to be the best way to deal with low/single-cell inputs

}
