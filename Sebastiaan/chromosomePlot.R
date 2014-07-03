# chomosomePlot.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of candidate gene data from Sebastiaan

chromosomePlot <- function(file = "BFMI861-S1vsALL_SNPs.txt"){
  chromosomes  <- c(1:19, "X", "Y", "M")
  
  chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
  snpOUT       <- read.table(file, sep="\t", header=TRUE)
  markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
  inconsistent <- read.table("Analysis/inconsistentSNPsAtlas.txt", sep="\t", header=TRUE)

  mlength <- max(chrInfo[,"Length"])
  plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="BFMI861-S1 versus the other BFMI",sub="Genotype errors: BFMI860-12, BFMI861-S2 and BFMI861-S1", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
  cnt <- 1
  aa <- apply(chrInfo,1,function(x){
    lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
    cnt <<- cnt + 1
  })
  aa <- apply(inconsistent, 1,function(x){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc, pch='|', col='black',cex=0.4)
  })
  aa <- apply(snpOUT, 1,function(x){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc+0.15, pch='▼', col='red',cex=0.5)
  })
  aa <- apply(markers, 1,function(x){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc-0.15, pch='▲', col='blue',cex=0.5)
  })
  axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
  axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
}

setwd("E:/Mouse/DiversityArray/")
chromosomePlot("Analysis/Diabetes/BFMI861-S1vsALL_SNPs.txt")
chromosomePlot("Analysis/Diabetes/BFMI861-S1vsBFMI861-S2_SNPs.txt")
chromosomePlot("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_SNPs.txt")

