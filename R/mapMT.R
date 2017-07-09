#' map the mitochondrial reads to the revised Cambridge reference sequence
#'
#' @param mtReads mitochondrial reads as a GAlignments, typically from getMT()
#' @param ...     other arguments to be passed to gsnap
#'
#' @import gmapR
#' @import rtracklayer
#' @import VariantTools
#'
#' @export
mapMT <- function(mtReads, ...) { 

  stop("mapMT() is not finished yet")

}

# Genome creation for gmapR:
#
# library(rtracklayer)
# rCRSfasta <- FastaFile(system.file("extdata", "NC_012920.1.fasta", 
#                                    package="ATACseeker", mustWork=TRUE))
# gmapGenomePath <- file.path("chrM_rCRS")
# gmapGenomeDirectory <- GmapGenomeDirectory(gmapGenomePath, create=TRUE)
# gmapGenomerCRS <- GmapGenome(rCRSfasta, gmapGenomeDirectory, "rCRS")

