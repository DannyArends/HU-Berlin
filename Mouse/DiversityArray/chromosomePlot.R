# chomosomePlot.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of candidate gene data from Sebastiaan

chromosomePlot <- function(file, main, plotMarkers = TRUE){
  chromosomes  <- as.character(c(1:19, "X", "Y", "M"))
  
  chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
  snpOUT       <- read.table(file, sep="\t", header=TRUE, colClasses=c("character"))
  markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
  inconsistent <- read.table("Analysis/inconsistentSNPsAtlas.txt", sep="\t", header=TRUE)

  mlength <- max(chrInfo[,"Length"])
  plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main=main, sub="Genotype errors: BFMI861-S2 and BFMI861-S1", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
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
    cat(x["Chr"],"->",yloc,"\n")
    points(x=xloc, y=yloc+0.15, pch='▼', col='red',cex=0.5)
  })
  if(plotMarkers){
    aa <- apply(markers, 1,function(x){
      yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
      points(x=xloc, y=yloc-0.15, pch='▲', col='blue',cex=0.5)
    })
  }
  axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
  axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
}

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")

png(file="Analysis/Figures/BFMI861-S1vsALL_SNPs.png",width=1200, height=800)
  chromosomePlot("Analysis/Diabetes/BFMI861-S1vsALL_SNPs.txt", "BFMI861-S1 versus the other BFMI")
dev.off()

png(file="Analysis/Figures/BFMI861-S1vsBFMI861-S2_SNPs.png",width=1200, height=800)
  chromosomePlot("Analysis/Diabetes/BFMI861-S1vsBFMI861-S2_SNPs.txt", "BFMI861-S1 versus BFMI861-S2")
dev.off()

png(file="Analysis/Figures/BFMI861-S1andBFMI860-12vsALL_SNPs.png",width=1200, height=800)
  chromosomePlot("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_SNPs.txt", "BFMI861-S1 equal to BFMI860-12 versus the other BFMI")
dev.off()

#Deike 1
png(file="Analysis/Figures/BFMI861-S2vs861S1n860-12n860-S2_SNPs.png",width=1200, height=800)
  chromosomePlot("Analysis/Diabetes/BFMI861-S2vs861S1n860-12n860-S2_SNPs.txt", "BFMI861-S2 different from BFMI861-S1,BFMI860-12, and BFMI860-S2", plotMarkers=FALSE)
dev.off()

#Deike 2
png(file="Analysis/Figures/BFMI861-S2nB6Nvs861S1n860-12n860-S2_SNPs.png",width=1200, height=800)
  chromosomePlot("Analysis/Diabetes/BFMI861-S2nB6Nvs861S1n860-12n860-S2_SNPs.txt", "BFMI861-S2 eq B6N different from BFMI861-S1,BFMI860-12, and BFMI860-S2", plotMarkers=FALSE)
dev.off()


#Deike 3
op <- par(mfrow=c(2,1))
png(file="Analysis/Figures/BFMI861-S1vs861S2n860-12n860-S2_SNPs.png",width=1200, height=800)
  chromosomePlot("Analysis/Diabetes/BFMI861-S1vs861S2n860-12n860-S2_SNPs.txt", "BFMI861-S1 different from BFMI861-S2,BFMI860-12, and BFMI860-S2", plotMarkers=FALSE)
dev.off()


