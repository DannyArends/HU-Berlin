# analysisOfMicrosatellitesInHorses.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Microsatellite Analysis

setwd("D:/Saria_Horses")
rawdata <- read.table("MicrosatellitesAnalyes.txt", sep='\t', row.names=1, header=TRUE, colClasses="character")

frequencies <- vector("list", ncol(rawdata))
names(frequencies) <- colnames(rawdata)
for(column in 1:ncol(rawdata)){
  genotypes <- unique(unlist(strsplit(unique(rawdata[,column]), "-")))
  genotypeData <- strsplit(rawdata[,column], "-")

  numbers <- NULL
  for(geno in genotypes){ numbers <- c(numbers, sum(unlist(lapply(genotypeData, function(x){x == geno})))) }

  freq <- (numbers / sum(numbers))
  names(freq) <- genotypes
  frequencies[[column]] <- freq
}

# Calculate the polymorphic information content, based on the genotype frequencies
PIC <- function(allelfreq){
  t1 <- 0
  for(Pi in allelfreq){ t1 = t1 + (Pi^2) }      # Sum of squared allele frequencies
  t2 <- 0
  for(i in 1:(length(allelfreq)-1)){
    for(j in (i+1):length(allelfreq)){
      Pi = allelfreq[i]; Pj = allelfreq[j]
      t2 = t2 + (Pi^2 * Pj^2)                   # Sum of squared between allele frequencies
    }
  }
  return(1 - t1 - 2 * t2)                       # PIC formula
}
