# differentialExpressionFigure.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
#
# Create a figure for the RNA sequencing data (pre-processed by MDC)

setwd("E:/Mouse/DNA/DiversityArray/")

chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
  
setwd("E:/Mouse/RNA/Sequencing/MPI_RPKM_ANALYSIS")

RPKM <- read.table("BFMI_RPKM_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")

mlength <- max(chrInfo[,"Length"])
plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="TestPlot", sub="Genotype errors: BFMI861-S2 and BFMI861-S1", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

#aa <- apply(RPKM, 1,function(x){
#  yloc <- match(as.character(x["chromosome_name"]), chromosomes); xloc <- as.numeric(x["start_position"])
#  ratio <- as.numeric(as.numeric(x["Ratio.1"]) < 1)
#  ttest <- as.numeric(x["tTest.1"])
#  col <- "white"
#  if(!is.na(ttest) && ttest < 0.1)  col <- rgb(.8, .8, .8)
#  if(!is.na(ttest) && ttest < 0.05) col <- rgb(.4, .4, .4)
#  if(!is.na(ttest) && ttest < 0.01) col <- rgb(0,0,0)
#  if(col != "white") points(x=xloc, y=yloc-0.15+(ratio*0.40), pch='|', col=col,cex=0.4)
#})

aa <- apply(RPKM, 1,function(x){
  yloc <- match(as.character(x["chromosome_name"]), chromosomes); xloc <- as.numeric(x["start_position"])
  #ttest <- as.numeric(x["tTest.1"])
  #if(!is.na(ttest) && ttest < 0.1){
    col <- "white"
    if(x["A.D_BFMI860.12xB6N"] == "B6N") col <- "gray"
    if(x["A.D_BFMI860.12xB6N"] == "BFMI") col <- "orange"
    if(x["A.D_BFMI860.12xB6N"] == "ADDITIVE") col <- "white"
    if(col != "white") points(x=xloc, y=yloc+0.2, pch='▲', col=col,cex=0.4)
    col <- "white"
    if(x["A.D_B6NxBFMI860.12"] == "B6N") col <- "gray"
    if(x["A.D_B6NxBFMI860.12"] == "BFMI") col <- "orange"
    if(x["A.D_B6NxBFMI860.12"] == "ADDITIVE") col <- "white"
    if(col != "white") points(x=xloc, y=yloc-0.08, pch='▼', col=col,cex=0.4)
  #}
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
