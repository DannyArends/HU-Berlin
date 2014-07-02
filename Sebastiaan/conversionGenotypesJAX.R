# Conversion of JAX Data to numerical genotypes
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#

setwd("E:/Mouse/DiversityArray/")

verbose <- FALSE
chromosomes <- c(1:19, "X", "Y", "MT")
rawchrdata <- NULL
for(chr in chromosomes){
  chrdata <- read.table(paste0("JAX/RAW/2009-11-29/Chr",chr,".txt"), sep="\t", header=TRUE, colClasses=c("character"))
  rawchrdata <- rbind(rawchrdata, chrdata)
}
genotypecolumns <- 14:ncol(rawchrdata)
nprobesOnChr <- nrow(rawchrdata)

cat("Loaded in", nprobesOnChr, "JAX genotypes on", length(genotypecolumns),"animals\n")

cat("JAX_ID", "\t", paste0(colnames(rawchrdata)[genotypecolumns], collapse="\t"), "\n", file="Analysis/MeasurementsJAX.txt", sep="")
aa <- lapply(1:nprobesOnChr, function(snp){
  A0    <- as.character(rawchrdata[snp, "Allele.A"])
  A2    <- as.character(rawchrdata[snp, "Allele.B"])
  Hetro <- "H"

  originalgenotypes <- as.character(rawchrdata[snp,genotypecolumns])
  correctedgenotypes <- rep(NA, length(genotypecolumns))

  is0 <- grep(A0, originalgenotypes)
  is2 <- grep(A2, originalgenotypes)
  isH <- grep(Hetro, originalgenotypes)

  wrongHin0 <- which(is0 %in% is2)
  wrongHin2 <- which(is2 %in% is0)
  isH <- c(isH, is0[which(is0 %in% is2)])
  if(length(wrongHin0) > 0) is0 <- is0[-wrongHin0]
  if(length(wrongHin2) > 0) is2 <- is2[-wrongHin2]

  if(length(is0) > 0) correctedgenotypes[is0] <- 0
  if(length(isH) > 0) correctedgenotypes[isH] <- 1
  if(length(is2) > 0) correctedgenotypes[is2] <- 2
  if(verbose) cat(paste0(A0,"/", A2,":"), originalgenotypes, "->", correctedgenotypes,"\n")
  if((snp %% 1000) == 0) cat("transformed ",snp,"/",nprobesOnChr,"\n", sep="")
  cat(rawchrdata[snp,1], "\t", paste0(correctedgenotypes, collapse="\t"), "\n", file="Analysis/MeasurementsJAX.txt", sep="", append=TRUE)
  return(NULL)
})

