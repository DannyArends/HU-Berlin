# analysisOfMicrosatellitesInHorses.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Microsatellite Analysis

setwd("D:/Saria_Horses")
rawdata <- read.table("MicrosatellitesAnalyes.txt", sep='\t', row.names=1, header=TRUE, colClasses="character")

rassen <- rawdata[,"Rassen"]
rawdata <- rawdata[,-c(1)]

difrassen <- unique(rassen)


calculateFrequencies <- function(microsatellites){
  frequencies <- vector("list", ncol(microsatellites))
  names(frequencies) <- colnames(microsatellites)
  for(column in 1:ncol(microsatellites)){
    genotypes <- unique(unlist(strsplit(unique(microsatellites[,column]), "-")))
    genotypeData <- strsplit(microsatellites[,column], "-")

    numbers <- NULL
    for(geno in genotypes){ numbers <- c(numbers, sum(unlist(lapply(genotypeData, function(x){x == geno})))) }

    freq <- (numbers / sum(numbers))
    names(freq) <- genotypes
    frequencies[[column]] <- freq
  }
  return(frequencies)
}

fAll <- calculateFrequencies(rawdata)
fH <- calculateFrequencies(rawdata[rassen=="H",])
fS <- calculateFrequencies(rawdata[rassen=="S",])
fK <- calculateFrequencies(rawdata[rassen=="K",])

numberOfGenotypes <- cbind(lapply(fAll,length),lapply(fH,length),lapply(fS,length),lapply(fK,length))
colnames(numberOfGenotypes) <- c("All", "H", "S", "K")


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
  return(as.numeric(1 - t1 - 2 * t2))                       # PIC formula
}

pAll <- lapply(fAll,PIC)
pH <- lapply(fH,PIC)
pS <- lapply(fS,PIC)
pK <- lapply(fK,PIC)


PICsPerRassen <- cbind(pAll,pH,pS,pK)


p <- NULL
for(x in seq(0,1, 0.01)){
 y <- 1-x
 p <- c(p, PIC(c(x,y)))
}
 
 
# HardyWeinberg:

datafor <- which(rawdata[,2] != "")
freq <- table(rawdata[datafor,2]) / length(datafor)

freqM = (0.08333333)^2 + 0.5 * 0.16666667
freqN = 1 - freqM

expectedMM = freqM^2
expectedMN = freqM * freqN
expectedNN = freqN^2