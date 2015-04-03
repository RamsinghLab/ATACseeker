## from Patrick Abouyoun's 2010 seminar
##
## e.g. for visualizing a specific enhancer:
##
## SMOenh <- GRanges('chr6', IRanges(29744702, 29751934), '+')
## seqinfo(SMOenh) <- seqinfo(Mus.musculus)[seqlevels(SMOenh)]
##
## plotCoverage(galp, SMOenh)
## plotCoverage(galp, SMOenh, stranded=TRUE)
## plotCoverage(get5primeCuts(galp), SMOenh, stranded=TRUE)
## plotCoverage(get5primeCuts(galp, shrink=FALSE), SMOenh, stranded=TRUE)
##
plotCoverage <- function(x, gr, col="violet", xlab="base", ylab="reads", 
                         main="Coverage", stranded=F, smoothen=0, kf='normal',
                         ...) {

  if(length(gr) > 1) {
    stop("This function only plots coverage for a single range...")
  } else {
    chr = as(seqnames(gr), 'character')
    if(smoothen > 0) {
      if(kf == 'exp') {
        wexp <- dexp(seq(0, 5, length.out=ceiling(smoothen/2)))
        wts <- round(c(rev(wexp), wexp)[1:smoothen], 5)
        wts <- wts/sum(wts)
      } else {
        wnorm <- dnorm(seq(0, 3, length.out=ceiling(smoothen/2)))
        wts <- round(c(rev(wnorm), wnorm)[1:smoothen], 5)
        wts <- wts/sum(wts)
      }
    }
    if(stranded == TRUE) { 
      pluscvg <- coverage(plusStrand(x), ...)[[chr]]
      minuscvg <- coverage(minusStrand(x), ...)[[chr]]
      revWin <- window(minuscvg, start(gr), end(gr))
      if(smoothen > 0) {
        revDepth = -1 * runwtsum(revWin, smoothen, wt=wts, endrule="constant")
      } else { 
        revDepth <- -1 * as(revWin, 'vector')
      }
    } else { 
      pluscvg <- coverage(x, ...)[[chr]]
      revDepth <- 0
    }
    fwdWin <- window(pluscvg, start(gr), end(gr))
    if(smoothen > 0) {
      fwdDepth = runwtsum(fwdWin, smoothen, wt=wts, endrule="constant")
    } else { 
      fwdDepth <- as(fwdWin, 'vector')
    }
  }

  x <- start(gr):end(gr)
  xlim <- c(start(gr), end(gr))
  ylim <- c(min(revDepth), max(fwdDepth))

  plot(x=start(gr), y=0, xlim=xlim, ylim=ylim, 
       xlab=xlab, ylab=ylab, main=main, type="n")
  polygon(c(start(gr), x, end(gr)), 
          c(0, as(fwdDepth, 'vector'), 0), 
          col=col)
  if(stranded == TRUE) {
    polygon(c(start(gr), x, end(gr)), 
            c(0, as(revDepth, 'vector'), 0), 
            col=col)
  }
}

## convenience function
plotStranded <- function(x, gr, ...) plotCoverage(x, gr, ..., stranded=TRUE)

