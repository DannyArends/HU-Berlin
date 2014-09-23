# DE_ASE_plots.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Sep, 2014

chromosomes  <- as.character(c(1:19, "X", "Y", "MT"))

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)

mlength <- max(chrInfo[,"Length"])

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
RPKM    <- read.csv("Analysis/BFMI_RPKM_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")
matBFMI <- read.csv("Imprinted_matBFMIsnps_5reads.txt", sep="\t", header=TRUE, colClasses="character")
matB6N  <- read.csv("Imprinted_matB6Nsnps_5reads.txt", sep="\t", header=TRUE, colClasses="character")

### maternal BFMI

plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="ASE & Dominant expression (maternal BFMI)", yaxt="n", xlab="Chromosome", ylab="Length (Mb)", xaxt="n")
abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

aa <- apply(matBFMI, 1,function(x){
  xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
  if(x["ImprintingScore"] > 0.35){
    col <- "gray"
    if(x["Origin"]=="BFMI") col <- "orange"
    points(x=xloc - 0.17, y=yloc, pch="-", col=col, cex=1.3)
  }
})

aa <- apply(RPKM, 1,function(x){
  xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
  col <- "white"
  if(x["A.D_BFMI860.12xB6N"] == "B6N") col <- "gray"
  if(x["A.D_BFMI860.12xB6N"] == "BFMI") col <- "orange"
  if(col != "white") points(x=xloc + 0.17, y=yloc, pch="-", col=col,cex=1.3)
})

cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1,lwd=3)
  cnt <<- cnt + 1
})

axis(1, chromosomes, at=c(1:nrow(chrInfo)), las=1)
axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
legend("topright", c("BFMI like", "B6N like"), fill=c("orange","gray"))

genesDominantExpressed <- RPKM[which(RPKM[,"A.D_BFMI860.12xB6N"] == "B6N" | RPKM[,"A.D_BFMI860.12xB6N"] == "BFMI"),"ensembl_gene_id"]
genesImprinted <- matBFMI[which(matBFMI[,"ImprintingScore"] > 0.35),"ensembl_gene_id"]
genesBoth <- genesDominantExpressed[which(genesDominantExpressed %in% genesImprinted)]
matBFMI[which(matBFMI[,"ensembl_gene_id"] %in% genesBoth),]

### maternal B6N

plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="ASE & Dominant expression (maternal BFMI)", yaxt="n", xlab="Chromosome", ylab="Length (Mb)", xaxt="n")
abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

aa <- apply(matB6N, 1,function(x){
  xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
  if(x["ImprintingScore"] > 0.35){
    col <- "gray"
    if(x["Origin"]=="BFMI") col <- "orange"
    points(x=xloc - 0.17, y=yloc, pch="-", col=col, cex=1.3)
  }
})

aa <- apply(RPKM, 1,function(x){
  xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
  col <- "white"
  if(x["A.D_B6NxBFMI860.12"] == "B6N") col <- "gray"
  if(x["A.D_B6NxBFMI860.12"] == "BFMI") col <- "orange"
  if(col != "white") points(x=xloc + 0.17, y=yloc, pch="-", col=col,cex=1.3)
})

cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1,lwd=3)
  cnt <<- cnt + 1
})

axis(1, chromosomes, at=c(1:nrow(chrInfo)), las=1)
axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
legend("topright", c("BFMI like", "B6N like"), fill=c("orange","gray"))

genesDominantExpressed <- RPKM[which(RPKM[,"A.D_B6NxBFMI860.12"] == "B6N" | RPKM[,"A.D_B6NxBFMI860.12"] == "BFMI"),"ensembl_gene_id"]
genesImprinted <- matB6N[which(matB6N[,"ImprintingScore"] > 0.35),"ensembl_gene_id"]
genesBoth <- genesDominantExpressed[which(genesDominantExpressed %in% genesImprinted)]
matB6N[which(matB6N[,"ensembl_gene_id"] %in% genesBoth),]