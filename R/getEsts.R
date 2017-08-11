#' Complexity estimation for ATACseq libraries.
#'
#' Estimate the complexity of a library or sample based on unique fragments
#' using Daley and Smith's implementation of Good-Toulmin rational function
#' approximation to solve the missing species problem.
#'
#' @param xx      The fragments or sample of fragments
#' @param withCI  Have preseq compute 95 percent confidence intervals for plots?
#' @param ...     Other arguments to pass on to preseq 
#' 
#' @return        A data frame with results
#' 
#' @export 
getEsts <- function(xx, withCI=FALSE, ...) {  
  if (class(xx) == "GRanges") xx <- getComplexity(xx)
  message("Estimating complexity (this can take a little while)...")
  if (!withCI) {
    res <- as.data.frame(preseqR.rfa.curve(xx, ...)$estimates)
    names(res) <- c("reads","frags")
  } else {
    res <- as.data.frame(preseqR.rfa.species.accum.curve(xx, ..., bootstrap=99))
    names(res) <- c("reads","frags","frags.lower","frags.upper")
  }
  return(res)
}
