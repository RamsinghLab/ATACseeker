#' what you would expect.  not terribly novel; maybe use esATAC's instead?
#' FIXME: generalize this to handle CTCF/TSS plots, etc. according to genome(x)
#'
#' @param   x         a GAlignmentPairs or GRanges to be scored w/featureScores
#' @param   anno      annotations, e.g. Homo.sapiens or Mus.musculus
#' @param   ordering  arg to pass to binPlots() 
#' @param   span      smoother span
#' @param   smoothing smoothing width
#' @param   ...       args to pass along to binPlots
#' 
#' @return  whatever binPlots returns
#' 
#' @import  Repitools
#'
#' @export
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
