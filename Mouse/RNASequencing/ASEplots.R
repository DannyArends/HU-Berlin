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
alldata <- read.csv("ReAnalysisSNPs/RPKM+ASE_min100reads_Threshold50.txt", sep="\t", header=TRUE, stringsAsFactor=FALSE, colClasses=c("character"))
mlength <- max(chrInfo[,"Length"])

png("ReAnalysisSNPs/ASE_SNPwise_Filtered100reads.png", width = 2000, height = 1000)
  op <- par(cex = 2.0, cex.main=2.0, cex.sub=1.3)
  plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="Allele specific expression", sub="paternal BFMI (left, BFMI♂ x B6N♀) versus maternal BFMI (right, B6N♂ x BFMI♀)", yaxt="n", xlab="", ylab="", xaxt="n")
  abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

  aa <- apply(alldata, 1, function(x){
    xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])

    col <- "white"
    if(x["matBFMI"]=="BFMI") col <- "gold1"
    if(x["matBFMI"]=="B6N") col <- "gray60"
    if(col != "white") points(x=xloc+0.2, y=yloc, pch="-", col=col, cex=1.8)
    
    col <- "white"
    if(x["matB6N"]=="BFMI") col <- "gold1"
    if(x["matB6N"]=="B6N") col <- "gray60"
    if(col != "white") points(x=xloc-0.2, y=yloc, pch="-", col=col, cex=1.8)
  })

  cnt <- 1
  aa <- apply(chrInfo,1,function(x){
    lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="black", lty=1, lwd=2)
    cnt <<- cnt + 1
  })

  axis(1,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1, cex.axis=1.5)
  axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=1.2, las=1)
  legend("topright", c("BFMI allele expressed", "B6N allele expressed"), fill=c("gold1","gray60"), cex=1.2)
dev.off()