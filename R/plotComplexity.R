#' plot results from preseq (complexity estimation around features)
#' TODO: use the ... to allow for comparison of multiple libraries
#'
#' @param   x       the thing whose complexity is to be plotted
#' @param   ...     other parameters to pass on to getEsts()
#' @param   ests    precomputed complexity estimates (speeds up multiple runs) 
#' @param   withCI  plot confidence intervals for the estimates?  (FALSE) 
#' @param   add     add to an existing plot? (FALSE) 
#' 
#' @return  invisibly, the complexity estimates corresponding to the plot
#' 
#' @seealso preseqR
#' 
#' @export
plotComplexity <- function(x, ..., ests=NULL, withCI=FALSE, add=FALSE) {

  if (!is.null(ests)) { 
    if (withCI & !all(c("frags.lower","frags.upper") %in% names(ests))) {
      stop("CIs requested, but provided estimates do not contain them!")
    }
    message("Using supplied estimates...")
    res <- ests
  } else if (is(x, "data.frame")) {
    res <- getEsts(x, ..., withCI=withCI)
  } else if (is(x, "GRanges") | is(x, "GAlignmentPairs")) {
    res <- getEsts(getComplexity(x), ..., withCI=withCI)
  } else {
    stop("Don't know what to do with a", class(x))
  }

  if (is.null(res)) {
    stop("Complexity estimate failed. You might try estimating from 5' cuts.")
  } else {
    message("Plotting complexity curve", 
           ifelse(withCI, "with confidence intervals...", "..."))
    if (withCI) { 
      # {{{
      epsilon <- 100
      if (add == TRUE) { 
        with(res, 
             points(x=reads, y=frags, pch=NA_integer_, ...))
      } else {
        with(res, 
             plot(x=reads, y=frags, 
                  main="Extrapolated library complexity (95% CI)", 
                  pch=NA_integer_, 
                  xlab="Read count", 
                  ylab="Unique fragments"))
      }
      with(res, 
           segments(reads - epsilon, frags.lower, 
                    reads + epsilon, frags.upper)) ## conf. int.
      points(x=res$reads, y=res$frags, pch=18) ## mean fit 
      lines(x=res$reads, y=res$frags, lwd=3, lty=1, col="red") ## mean fit 
      lines(x=res$reads, y=res$frags.lower, lwd=1, lty=3, col="grey50")
      lines(x=res$reads, y=res$frags.upper, lwd=1, lty=3, col="grey50")
      # }}}
    } else { 
      # {{{
      if (add == TRUE) { 
        with(res, 
             points(x=reads, y=frags, pch=NA_integer_, ...))
      } else {
        with(res,
             plot(x=reads, y=frags, 
                  col="red", pch=18,
                  main="Extrapolated library complexity",
                  xlab="Read count", 
                  ylab="Unique fragments"))
      }
      with(res,
           lines(x=reads, y=frags, lwd=3, lty=3, col="grey50"))
      # }}}
    }
  }
  invisible(res)

}

