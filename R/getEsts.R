#' Complexity estimation for ATACseq libraries.
#'
#' Estimate the complexity of a library or sample based on unique fragments
#' using Daley and Smith's implementation of Good-Toulmin rational function
#' approximation to solve the missing species problem.
#' 
#' Original functions in preseqR v2.0.1.1 for this were:
#' preseqR.rfa.curve and preseqR.rfa.species.accum.curve
#' 
#' The new functions as of the version 3.1.2 are:
#' 
#' ds.mincount == preseqR.rfa.curve
#' ds.mincount.bootstrap == preseqR.rfa.species.accum.curve
#' 
#' The new functions return generators that can have data passed to them
#' instead of returning a data frame as in version 2.0.1.1.
#' 
#' It is worth noting that these functions return *non-bootstrapped* 
#' fragment estimates but uses the bootstrapped variance estimates to 
#' calculate CI. This is a function of the latest implementation of preseqR.
#'
#' @param xx      The fragments or sample of fragments
#' @param withCI  Have preseq compute 95 percent confidence intervals for plots?
#' @param ...     Other arguments to pass on to preseqR 
#' 
#' @return        A data frame with results
#' 
#' @import preseqR
#' 
#' @export 
getEsts <- function(xx, withCI=FALSE, ...) {  
  if (is(xx, "GRanges") | is(xx, "GAlignmentPairs")) xx <- getComplexity(xx)
  message("Estimating complexity (this can take a little while)...")
  if (!withCI) {
    res.fun <- ds.mincount(xx, mt = 100, ...)
    message("Done generating complexity estimator...")
    message("Building library complexity curves...")
    res <- data.frame(reads = as.double(xx[, 1] %*% xx[, 2]) * 1:100,
                      frags = sapply(1:100, res.fun$FUN))
  } else {
    res.fun <- ds.mincount.bootstrap(xx, mt = 100, times = 99, ...)
    message("Done generating complexity estimator...")
    message("Building library complexity curves...")
    message("Calculating 95% CI...")
    res.ci <- sapply(1:100, function(x) {
      exp(qnorm((1 + 0.95) / 2.0) * sqrt(log(1.0 + res.fun$var(x) / (res.fun$FUN.nobootstrap$FUN(x)^2))))
    })
    res <- data.frame(reads = as.double(xx[, 1] %*% xx[, 2]) * 1:100,
                      frags = sapply(1:100, res.fun$FUN.nobootstrap$FUN))
    res$frags.lower = res$frags / res.ci
    res$frags.upper = res$frags * res.ci
  }
  return(res)
}
