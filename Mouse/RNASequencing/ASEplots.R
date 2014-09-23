# ASEplots.R - Analyse the SNPs and indels called by the GenomeAnalysisToolKit
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sep, 2014
# first written Sep, 2014

# ASE plots
chromosomes  <- as.character(c(1:19, "X", "Y", "MT"))

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
matBFMI <- read.csv("Imprinted_matBFMIsnps_5reads.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)
matB6N <- read.csv( "Imprinted_matB6Nsnps_5reads.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE)

mlength <- max(chrInfo[,"Length"])

plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="Origin of imprinted genes", yaxt="n", xlab="Chromosome", ylab="Length (Mb)", xaxt="n")
abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

aa <- apply(matBFMI, 1,function(x){
  xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
  if(x["ImprintingScore"] > 0.35){
    col <- "gray"
    if(x["Origin"]=="BFMI") col <- "orange"
    points(x=xloc+0.2, y=yloc, pch="-", col=col, cex=1.3)
  }
})
  
aa <- apply(matB6N, 1,function(x){
  xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
  if(x["ImprintingScore"] > 0.35){
    col <- "gray"
    if(x["Origin"]=="BFMI") col <- "orange"
    points(x=xloc-0.2, y=yloc, pch="-", col=col, cex=1.3)
  }
})

cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1, lwd=2)
  cnt <<- cnt + 1
})

axis(1,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
legend("topright", c("maternal BFMI", "maternal B6N"), fill=c("orange","gray"))
