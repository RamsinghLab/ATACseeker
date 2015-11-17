#' Plot fragment size density for a set of samples.
#' 
#' @param frag.dists    the computed densities of each fragment distribution
#' @param label         the group label to be used in the title
#'
#' @seealso fragmentLengths
#' 
#' @export
fragPlot <- function(frag.dists, label) {
  dens.max <- max(do.call(c, lapply(frag.dists, `[[`, "y")))
  plot(frag.dists[[1]], lwd=3, ylim=c(0, dens.max * 1.1),
       xlab="Fragment sizes (base pairs)", ylab="Density", 
       main=paste(label, "sample fragment sizes"), col=1)
  for (i in seq_along(frag.dists)[-1]) {
    lines(frag.dists[[i]], col=i, lwd=3)
  }
  getMode <- function(dens) dens$x[which.max(dens$y)]
  frag.mode <- median(round(do.call(c, lapply(frag.dists, getMode))))
  abline(v=frag.mode, col="gray80", lwd=3, lty=3)
  legend("topright", legend=names(frag.dists), col=seq(length(frag.dists)), 
         lwd=1, bty="n")
}
