##
## for getComplexity and preseqR 
##
getFeatureCounts <- function(x) {

  if (!is(x, 'GRanges')) x <- as(x, 'GRanges')
  ends <- paste(start(x), end(x), sep=':')
  freqs <- data.frame(seen=as(table(ends), 'matrix')[,1])
  freqs[ rev(order(freqs$seen)), , drop=FALSE] ## descending order 

}

##
## formatted for preseqR 
##
getComplexity <- function(x) {

  if (is(x, 'GAlignmentPairs')) x <- getFeatureCounts(as(x, 'GRanges'))
  if (!'seen' %in% names(x)) x <- getFeatureCounts(x)
  freq <- seq_len(max(range(x$seen)))
  breaks <- c(0, freq)
  d <- data.frame(freq=freq, species=hist(x$seen, breaks=breaks, plot=F)$counts)
#  d <- d[ which(d[,2] != 0), ] 
  return(data.matrix(d))

}

## 
## preseq applied to the above
##
## FIXME: use the ... to allow for comparison of multiple libraries
##
plotComplexity <- function(x, ..., withCI=FALSE) {

  ## require(preseqR)

  getEsts <- function(xx, withCI=withCI, ...) {  
    message("Nota bene: this can take a little while.")
    if (!withCI) {
      res <- as.data.frame(preseqR.rfa.curve(xx, ...)$estimates)
      names(res) <- c('reads','frags')
    } else {
      res <- as.data.frame(preseqR.rfa.species.accum.curve(xx, ...))
      names(res) <- c('reads','frags','frags.lower','frags.upper')
    }
    return(res)
  } 

  if (is(x, 'data.frame')) {
    res <- getEsts(x, withCI=withCI, ...)
  } else if (is(x, 'GRanges') | is(x, 'GAlignmentPairs')) {
    res <- getEsts(getComplexity(x), withCI=withCI, ...)
  } else {
    stop("Don't know what to do with a", class(x))
  }

  if (is.null(res)) {
    stop("Complexity estimate failed. You might try estimating from 5' cuts.")
  } else {
    message("Plotting complexity curve", 
           ifelse(withCI, "with confidence intervals", ""))
    if (withCI) { 
      epsilon <- 100
      plot(x=res$reads, y=res$frags, 
           main="Extrapolated library complexity (95% CI)", 
           pch=NA_integer_, xlab="Unique fragments", ylab="Fragment count")
      segments(res$reads-epsilon, res$frags.lower, 
               res$reads+epsilon, res$frags.upper) ## conf. int.
      points(x=res$reads, y=res$frags, pch=18) ## mean fit 
      lines(x=res$reads, y=res$frags, lwd=3, lty=1, col="red") ## mean fit 
      lines(x=res$reads, y=res$frags.lower, lwd=1, lty=3, col="grey50")
      lines(x=res$reads, y=res$frags.upper, lwd=1, lty=3, col="grey50")
    } else { 
      plot(x=res$reads, y=res$frags, col='red', pch=18,
           main="Extrapolated library complexity",
           xlab="Unique fragments", ylab="Fragment count")
      lines(x=res$reads, y=res$frags, lwd=3, lty=3, col="grey50")
    }
  }
  invisible(res)

}
