# RNA Seq - Expression data analysis, comparison back to MPI data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2014
# first written Aug, 2014

## Comparison to MPI (Broken because of weird names)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

RPKM_MPI <- read.csv("MPI_RPKM_ANALYSIS/2014-07-04_BFMI_RPKM.txt", sep="\t", header=TRUE, row.names=1)       # RNA-Seq summary data from MDC
cat("MPI called expressions for", nrow(RPKM_MPI), "genes\n")

corenames <- as.character(sampleIDs[match(as.numeric(unlist(lapply(strsplit(colnames(RPKM),"_"),"[",1))), sampleIDs$Lib_id), "core_name"])
corenames <- unlist(lapply(strsplit(corenames,"_"),"[",1))

shortnames <- unlist(lapply(strsplit(colnames(RPKM),"_"),"[",1))

RPKM_MPI <- RPKM_MPI[, unlist(lapply(corenames, grep, colnames(RPKM_MPI)))]

RPKM_MPI <- RPKM_MPI[which(rownames(RPKM_MPI) %in% rownames(RPKM)),]
RPKM_GAB <- RPKM[which(rownames(RPKM) %in% rownames(RPKM_MPI)),]

png("comparison_RPKM_GAB_MPI.png")
  op <- par(mfrow=c(1,2))
  plot(RPKM_GAB[,1], RPKM_MPI[,1], main= colnames(RPKM_GAB)[1], log = "xy")
  plot(RPKM_GAB[,2], RPKM_MPI[,2], main= colnames(RPKM_GAB)[2], log = "xy")
dev.off()
