## e.g.
##
## data(cytobands_hg19)
## del7q <- byBand(cytobands_hg19)$chr7q22.1
## galp <- atacPairedEnd('chr7.q10.CD49f_r1.hg19.bam', which=del7q)
##
## ...or...
##
## library(Mus.musculus)
## SMO <- transcriptsBy(Mus.musculus, 'gene')[['319757']]
## galp <- atacPairedEnd('A1.SMO.bam', which=SMO)
##
## ...or...
##
## library(Mus.musculus)
## chr6 <- GRanges('chr6', IRanges(1, seqlengths(Mus.musculus)['chr6']), '*')
## galp <- atacPairedEnd("A2.mm10.unique.bam",
##                       bamParams=properPairedEndAtacFilters(which=chr6))
##
atacPairedEnd <- function(bam, genome=c('hg19','mm10'), bamParams=NULL, which=NULL, ...) {

  getStdChromGRanges <- function(x) {
    ## ONLY works if chromosomes are properly ordered as in OrganismDbi
    as(seqinfo(x), 'GRanges')[ 1:(which(seqlevels(x) == 'chrM') - 1) ] 
  }

  if(is.null(bamParams)) {
    if (is.null(which)) {
      genome <- match.arg(genome)
      which <- switch(genome,
                      hg19 = getStdChromGRanges(Homo.sapiens),
                      mm10 = getStdChromGRanges(Mus.musculus))
    } 
    bamParams <- properPairedEndAtacFilters(which=which, ...)
  }

  readGAlignmentPairsFromBam(bam, param=bamParams)

}

## not exported; used for preseq estimation
properPairedEndAtacFilters <- function(which, ...) {

  ScanBamParam(what=c('rname','strand','pos','isize','mapq'),
               flag=scanBamFlag(isProperPair=TRUE,
                                isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

## not exported; used for BAM filtering, also needs which(!chrM) filter (duh?)
uniquePairedEndAtacFilters <- function(which, ...) {

  ScanBamParam(what=c('rname','strand','pos','isize','mapq'),
               flag=scanBamFlag(isProperPair=TRUE,
                                isDuplicate=FALSE,
                                isNotPrimaryRead=FALSE,
                                isNotPassingQualityControls=FALSE), 
               which=which, ...)
}

## for recovering spike-ins
spikeInFilters <- function(...) {

  ## grab the unmapped sequences for brute-force alignment to phiX spike-ins
  ScanBamParam(what=c('seq'), 
               flag=scanBamFlag(isUnmappedQuery=TRUE, 
                                isNotPassingQualityControls=FALSE, 
                                ...))
}
