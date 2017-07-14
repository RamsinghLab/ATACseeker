# how the mtGenes GRanges were generated:
if (FALSE) {   
  mtBp <- 16569
  mtChr <- "chrM"
  mtRanges.hg19 <- list("Control-Region"=c(1,576),
                        "tRNA"=c(577,647),
                        "rRNA"=c(648,1601), 
                        "tRNA"=c(1602,1670),
                        "rRNA"=c(1671,3229),
                        "tRNA"=c(3230,3304), 
                        "MT-ND1"=c(3307,4262),
                        "tRNA"=c(4263,4400),
                        "tRNA"=c(4402,4469),
                        "MT-ND2"=c(4470,5511),
                        "tRNA"=c(5512,5579),
                        "tRNA"=c(5587,5655),
                        "tRNA"=c(5657,5891),
                        "MT-CO1"=c(5904,7445),
                        "tRNA"=c(7446,7514),
                        "tRNA"=c(7518,7585),
                        "MT-CO2"=c(7586,8269),
                        "tRNA"=c(8295,8364),
                        "MT-ATP8"=c(8366,8572),
                        "MT-ATP6"=c(8573,9207),
                        "MT-CO3"=c(9208,9990),
                        "tRNA"=c(9991,10058),
                        "MT-ND3"=c(10059,10404),
                        "tRNA"=c(10405,10469),
                        "MT-ND4L"=c(10470,10766),
                        "MT-ND4"=c(10767,12137),
                        "tRNA"=c(12138,12336),
                        "MT-ND5"=c(12337,14148),
                        "MT-ND6"=c(14149,14673),
                        "tRNA"=c(14674,14742),
                        "MT-CYB"=c(14747,15887),
                        "tRNA"=c(15888,15953),
                        "tRNA"=c(15956,16023),
                        "Control-Region"=c(16024,mtBp))
  mtGenes.hg19 <- suppressWarnings(data.frame(do.call(rbind, mtRanges.hg19)))
  names(mtGenes.hg19) <- c("start","end")
  mtGenes.hg19$chrom <- rep(mtChr, nrow(mtGenes.hg19))
  mtGenes.hg19 <- makeGRangesFromDataFrame(mtGenes.hg19) 
  mtGenes.hg19$name <- names(mtRanges.hg19)

  library(rtracklayer)
  hg19ToHg38 <- system.file("extdata", "hg19ToHg38.over.chain", 
                            package="ATACseeker", mustWork=TRUE)
  chain <- import.chain(hg19ToHg38)
  mend <- function(x) {
    start(x) <- min(start(x))
    end(x) <- max(end(x))
    reduce(x)
  }  
  mtGenes.hg38 <- unlist(GRangesList(lapply(liftOver(mtGenes.hg19,chain),mend)))
  mtGenes.hg38$name <- mtGenes.hg19$name
}
