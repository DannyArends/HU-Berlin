# InterestingRegions.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of candidate gene data from Sebastiaan

setwd("E:/Mouse/DNA/DiversityArray/")

snpOUT  <- read.table("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

snpInRegion <- 35
maxDistance <- 4000000

regions <- NULL
for(chr in chromosomes){
  onchr  <- snpOUT[snpOUT[,"Chr"] == chr,]
  if(nrow(onchr) > snpInRegion){
    sorted <- sort(as.numeric(onchr[,"Location"]), index.return=TRUE)
    onchr  <- onchr[sorted$ix,]
    sloc   <- 0
    snpcnt <- 0
    for(snp in 1:(nrow(onchr)-1)){
      cloc <- onchr[snp, "Location"]
      nloc <- onchr[snp+1, "Location"]
      if(as.numeric(nloc) - as.numeric(cloc) < maxDistance){
        if(snpcnt == 0) sloc <- cloc
        snpcnt <- snpcnt + 1
      }else{
        if(snpcnt > snpInRegion){
          regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
          cat("New region:", chr, sloc, cloc, snpcnt,"\n")
        }
        snpcnt <- 0
      }
    }
    if(snpcnt > snpInRegion){
      regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
      cat("New region:",chr, sloc, cloc, snpcnt,"\n")
    }
  }
}

colnames(regions) <- c("Chr","Start","End","SNPs")