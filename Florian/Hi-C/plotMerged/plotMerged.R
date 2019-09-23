#
# Plot results from a Hi-C experiment
# (c) Danny Arends (HU-Berlin) Sept - 2019
#

library(h5)
f <- h5file("C:/Users/Arends/Downloads/ENCFF971GRZ.h5")

setwd("D:/")

digested <- read.csv("HindIII_Digested.bed", sep="\t", header=FALSE)
options(scipen=999)

# Fix the chromosome name differences, and write out the genome_bins.txt file for binning
bins <- f["bin_positions"][]
ourbins <- f["bin_positions"][]
chromosomes <- as.character(unique(digested[,1]))
x <- 1
for(chr in unique(bins[,1])){
  ourbins[which(as.character(bins[,1]) == chr),1] <- chromosomes[x]
  x <- x+1
}
write.table(ourbins, "genome_bins.txt",sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

# Use the genome_bins.txt, and perform the binning on the server, the created file can be loaded into R
interactions <- read.table("binned.alignments.txt", sep = "\t")

# Plot comparing the first 100 bins on chromosome 1
op <- par(mfrow=c(2,1))
image(f["interactions"][1:100, 1:100], breaks = c(0, 10, 100, 10000), col=c("white", "gray", "black"))
box()
image(as.matrix(interactions[1:100, 1:100]), breaks = c(0, 10, 100, 10000), col=c("white", "gray", "black"))
box()

# Correlation between the two analysis paths
allC <- c()
for(x in 1:ncol(interactions)){
  allC <- c(allC, cor(interactions[,x], f["interactions"][,x], use = "pair"))
}

h5close(f)
