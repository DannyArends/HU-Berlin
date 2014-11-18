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
allData <- read.csv("ASE_10reads_combined.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE, colClasses=c("character"))
matBFMI <- allData[which(allData[,"cross"]=="matBFMI"),]
matB6N  <- allData[which(allData[,"cross"]=="matB6N"),]

mlength <- max(chrInfo[,"Length"])

plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="Allele specific expression", sub="maternal B6N (left, B6N♀ x BFMI♂) versus maternal BFMI (right, BFMI♀ x B6N♂)", yaxt="n", xlab="", ylab="", xaxt="n")
abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

aa <- apply(matBFMI, 1, function(x){
  xloc <- match(as.character(x["Chr"]), chromosomes); yloc <- as.numeric(x["Loc"])

  col <- "gray60"
  if(x["Origin1"]=="BFMI") col <- "gold1"
  points(x=xloc+0.2, y=yloc, pch="-", col=col, cex=2.9)

})
  
aa <- apply(matB6N, 1,function(x){
  xloc <- match(as.character(x["Chr"]), chromosomes); yloc <- as.numeric(x["Loc"])

    col <- "gray60"
    if(x["Origin1"]=="BFMI") col <- "gold1"
    points(x=xloc-0.2, y=yloc, pch="-", col=col, cex=2.9)

})

cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="black", lty=1, lwd=2)
  cnt <<- cnt + 1
})

axis(1,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1, cex.axis=1.5)
axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=1.2, las=1)
legend("topright", c("BFMI Allele expressed", "B6N Allele expressed"), fill=c("gold1","gray60"), cex=1.2)
