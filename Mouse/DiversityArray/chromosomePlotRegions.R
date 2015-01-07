# chomosomePlotRegions.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of candidate gene data from Sebastiaan

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
snpOUT       <- read.table("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
regions      <- read.table("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_MarkersInRegions.txt", sep="\t", header=TRUE, colClasses=c("character"))
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)

selectedmarkers <- unlist(lapply(regions[,"markersInRegions"], strsplit,", "))

mlength <- max(chrInfo[,"Length"])

png(file="Analysis/Figures/BFMI861-S1andBFMI860-12vsALL_SNPs_Regions.png",width=1200, height=800)

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="BFMI861-S1 equal to BFMI860-12 versus the other BFMI", sub="Genotype errors: BFMI861-S2 and BFMI861-S1", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(snpOUT, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  #cat(x["Chr"],"->",yloc,"\n")
  points(x=xloc, y=yloc+0.25, pch='▼', col='red',cex=0.5)
})

aa <- apply(markers, 1,function(x){
  if(x["markerID"] %in% selectedmarkers){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc-0.15, pch='▲', col='blue',cex=0.5)
  }
})

aa <- apply(regions, 1, function(x){
  yloc <- match(x["Chr"], chromosomes); xlocS <- x["Start"]; xlocE <- x["End"]
  lines(c(as.numeric(xlocS)-2500000,as.numeric(xlocE)+2500000), c(yloc, yloc), type="l", col="black", lty=1,lwd=5)
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)

dev.off()

snpOUT861S1  <- read.table("Analysis/Diabetes/BFMI861-S1vsALL_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
regions861S1 <- read.table("Analysis/Diabetes/BFMI861-S1vsALL_MarkersInRegions.txt", sep="\t", header=TRUE, colClasses=c("character"))

selectedmarkers861S1 <- unlist(lapply(regions861S1[,"markersInRegions"], strsplit,", "))

png(file="Analysis/Figures/RegionComparison_BFMI861-S1andBFMI860-12_vs_BFMI861-S1.png",width=1200, height=800)

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="Comparison of regions", sub="Genotype errors: BFMI861-S2 and BFMI861-S1", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(regions, 1, function(x){
  yloc <- match(x["Chr"], chromosomes); xlocS <- x["Start"]; xlocE <- x["End"]
  lines(c(as.numeric(xlocS)-2500000,as.numeric(xlocE)+2500000), c(yloc, yloc)+0.1, type="l", col="red", lty=1,lwd=5)
})

aa <- apply(regions861S1, 1, function(x){
  yloc <- match(x["Chr"], chromosomes); xlocS <- x["Start"]; xlocE <- x["End"]
  lines(c(as.numeric(xlocS)-2500000,as.numeric(xlocE)+2500000), c(yloc, yloc)-0.1, type="l", col="blue", lty=1,lwd=5)
})

aa <- apply(snpOUT, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  #cat(x["Chr"],"->",yloc,"\n")
  points(x=xloc, y=yloc+0.15, pch='|', col='yellow',cex=0.4)
})

aa <- apply(snpOUT861S1, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  #cat(x["Chr"],"->",yloc,"\n")
  points(x=xloc, y=yloc-0.05, pch='|', col='green',cex=0.4)
})

aa <- apply(markers, 1,function(x){
  if(x["markerID"] %in% selectedmarkers){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc-0.15, pch='▲', col='red',cex=0.7)
  }
  if(x["markerID"] %in% selectedmarkers861S1){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc+0.15, pch='▼', col='blue',cex=0.7)  
  }
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
legend("topright",c("BFMI861-S1 == BFMI860-12", "BFMI861-S1", "Marker (red region)", "Marker (blue region)","SNP (red region)", "SNP (blue region)"), lty=c(1,1,0,0,0,0), lwd=5, col=c("red","blue","red","blue","yellow","green"), pch=c("","","▲","▼","|","|"))

dev.off()




png(file="Analysis/Figures/BFMI861-S1vsALL_SNPs_Regions.png",width=1200, height=800)

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="BFMI861-S1 versus the other BFMI", sub="Genotype errors: BFMI861-S2 and BFMI861-S1", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(snpOUT861S1, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  #cat(x["Chr"],"->",yloc,"\n")
  points(x=xloc, y=yloc+0.25, pch='▼', col='red',cex=0.5)
})

aa <- apply(markers, 1,function(x){
  if(x["markerID"] %in% selectedmarkers861S1){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc-0.15, pch='▲', col='blue',cex=0.5)
  }
})

aa <- apply(regions861S1, 1, function(x){
  yloc <- match(x["Chr"], chromosomes); xlocS <- x["Start"]; xlocE <- x["End"]
  lines(c(as.numeric(xlocS)-2500000,as.numeric(xlocE)+2500000), c(yloc, yloc), type="l", col="black", lty=1,lwd=5)
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)

dev.off()

