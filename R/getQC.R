#' obtain QC information used by csaw from paired-end BAM files 
#' 
#' @param name    the name of the sample
#' @param pe.bam  the name of the BAM file
#'
#' @return a list with components "qc" (metrics) and "frag.dist" (density)
#'
#' @export
getQC <- function(name, pe.bam) {
  out <- getPESizes(pe.bam)
  too.large <- sum(out$sizes > 400)
  qc <- c(sample=name, out$diagnostics, too.large=too.large)
  frag.dist <- density(out$sizes[out$sizes <= 800])
  results <- list(qc=qc, frag.dist=frag.dist)
  return(results)
}
