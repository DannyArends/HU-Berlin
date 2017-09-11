# locations on the LWTL01 genome
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")

locLWLT01 <- read.table("FilteredLocationLWLT01.txt", header=TRUE, colClasses=c("character","character","numeric","numeric","character", "character","character"), row.names=1)
locLWLT01 <- cbind(locLWLT01, snpLoc = (locLWLT01[,"Start"] + locLWLT01[,"Stop"]) / 2)
locLWLT01 <- locLWLT01[with(locLWLT01, order(snpLoc)), ]

chrs <- as.character(1:29)

locLWLT01ordered <- NULL
for(chr in chrs){
  locLWLT01ordered <- rbind(locLWLT01ordered, locLWLT01[which(locLWLT01[, "chrN"] == chr),])
}

snpdataSiham <- read.table("filtered_snps.txt",sep="\t", check.names=FALSE)

locLWLT01ordered <- locLWLT01ordered[which(rownames(locLWLT01ordered) %in% rownames(snpdataSiham)),]

fullSNPdata <- cbind(locLWLT01ordered, snpdataSiham[rownames(locLWLT01ordered),])
fullSNPdata <- cbind(SNPname = rownames(fullSNPdata), fullSNPdata)


write.table(fullSNPdata, file="SNPnGeno_LWTL01.txt",sep="\t", row.names=FALSE)

setwd("D:/Edrive/Goat/DNA/Ibex")

ibex <- read.table("Goat_Nov2016_FinalReport.txt", sep = "\t",skip = 9, header=TRUE, row.names=1, check.names=FALSE)
samples <- read.table("Goat_Nov2016_Samples.txt", sep = "\t", header=TRUE)
colnames(ibex) <- samples[, "ID"]

fullIbexdata <- cbind(locLWLT01ordered, ibex[rownames(locLWLT01ordered),])
fullIbexdata <- cbind(SNPname = rownames(fullIbexdata), fullIbexdata)

write.table(fullIbexdata, file="SNPnGenoIBEX_LWTL01.txt",sep="\t", row.names=FALSE)