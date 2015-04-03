##
## e.g.
##
## galp <- data(humanGalp) 
## library(Homo.sapiens)
## anno <- promoters(TxDb.Hsapiens.UCSC.hg19.knownGene)
## tssQcPlot(galp, anno)
##
## ...or...
##
## galp <- data(mouseGalp) 
## library(Mus.musculus)
## data(txsByFpkmInMouseLskCells)
## anno <- flank(txsByFpkmInMouseLskCells, 200)
## anno <- anno[ which(anno$FPKMinLSK > 0) ] 
## tssQcPlot(galp, anno, mcols(anno))
## 
## FIXME: generalize this to handle CTCF/TSS plots, etc. according to genome(x)
## 
## data(CTCFsites.hg19)
## data(CTCFsites.mm10)
## 
tssQcPlot <- function(x, anno, ordering, span=5000, smoothing=50, ...) {

  strip <- function(x) { 
    mcols(x) <- NULL
    return(x)
  }
  
  ## require(Repitools)
  if (is(x, 'GAlignmentPairs')) {
    covs <- featureScores(as(c(strip(right(x)), strip(left(x))), 'GRanges'),
                          anno, up=span, down=span, freq=100, s.width=smoothing)
  } else if (is(x, 'GRanges')) {                          
    covs <- featureScores(strip(x), 
                          anno, up=span, down=span, freq=100, s.width=smoothing)
  }
  names(covs) <- 'ATACseq'
  binPlots(covs, ordering=ordering, ord.label="by quantile", n.bins=5, ...)

}
