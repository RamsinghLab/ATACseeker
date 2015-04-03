fragmentLengths <- function(galp, plotAs=c('hist','density','none'), ...){

  stopifnot(all(isProperPair(galp)))

  plotAs <- match.arg(plotAs, c('hist','density','none'))
  fraglens <- width(granges(galp))

  if( plotAs == 'hist' ) {  
    hist(fraglens, 
         main='Fragment size distribution', 
         ylab='paired-end reads of this size', 
         xlab='fragment length', 
         xlim=c(0, max(fraglens)),
         col='violet')
  } else if( plotAs == 'density' ) {
    d <- density(fraglens, ...)
    plot(d, zero.line=TRUE, 
         main='Fragment size distribution', 
         ylab='proportion of pairs with this fragment length', 
         xlab='fragment length', 
         xlim=c(0, max(fraglens)+100),
         type='s')
    polygon(d, col='violet', border='black')
  }

  invisible(fraglens)

}

