# RNA Seq - Expression data analysis, comparison back to MPI data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2014
# first written Aug, 2014

## Comparison to MPI (Broken because of weird names)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

sampleIDs <- read.table("FASTQ/sampledescription.txt",sep="\t", header=TRUE)
RPKM      <- read.csv("Analysis/RPKMnorm.txt", sep="\t", check.names=FALSE, header=TRUE, row.names=1)                             # RNA-Seq summary data
cat("Called expressions for", nrow(RPKM), "genes\n")                                                                              # Called expressions for 41388 genes
RPKM_MPI  <- read.csv("MPI_RPKM_ANALYSIS/2014-07-04_BFMI_RPKM.txt", check.names=FALSE, sep="\t", header=TRUE, row.names=1)        # RNA-Seq summary data from MDC
cat("MPI called expressions for", nrow(RPKM_MPI), "genes\n")                                                                      # MPI called expressions for 23330 genes

libraryIDs  <- as.numeric(unlist(lapply(strsplit(colnames(RPKM),"_"),"[",1)))
RPKM        <- RPKM[,which(libraryIDs %in% sampleIDs$Lib_id)]
libraryIDs  <- libraryIDs[which(libraryIDs %in% sampleIDs$Lib_id)]

orderInSampleID <- match(libraryIDs, sampleIDs$Lib_id)
corenames       <- as.character(sampleIDs[orderInSampleID, "core_name"])
RPKM            <- RPKM[,orderInSampleID]
corenames       <- unlist(lapply(lapply(strsplit(corenames,"_"),"[",c(1:2)), function(x){ paste0(x[1],"_",x[2]) }))

colnames(RPKM)      <- corenames
colnames(RPKM_MPI)  <- gsub("RPKM_","",colnames(RPKM_MPI))

RPKM_MPI <- RPKM_MPI[, unlist(lapply(corenames, grep, colnames(RPKM_MPI)))]

RPKM_MPI <- RPKM_MPI[which(rownames(RPKM_MPI) %in% rownames(RPKM)),]
RPKM_GAB <- RPKM[which(rownames(RPKM) %in% rownames(RPKM_MPI)),]

png("comparison_RPKM_GAB_MPI.png")
  op <- par(mfrow=c(1,2))
  plot(RPKM_GAB[,1], RPKM_MPI[,1], main= colnames(RPKM_GAB)[1], log = "xy")
  plot(RPKM_GAB[,2], RPKM_MPI[,2], main= colnames(RPKM_GAB)[2], log = "xy")
dev.off()
