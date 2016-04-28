# Preprocessing of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016

setwd("E:/Goat/DNA/Annotation")

arrayinfo <- read.csv("IGGC CapriBatch1.txt", sep="\t", skip = 29)
arrayinfo[1:10,]
posUncertain <- as.character(arrayinfo[which(duplicated(arrayinfo[,"loc_snp_id"])),2])

arrayinfo <- cbind(arrayinfo, reference = NA)

for(chr in c(1:29, "X")){
  dnasequence <- readLines(gzfile(paste0("chr",chr,".fna.gz")))
  dnasequence <- paste0(dnasequence[-1], collapse="")
  dnasequence <- strsplit(dnasequence, "")[[1]]
  arrayinfo[arrayinfo[,"chr"] == chr, "reference"] <- dnasequence[arrayinfo[arrayinfo[,"chr"] == chr,"chr_pos"]]
}

setwd("E:/Goat/DNA/Siham RawData")

ILMNannot <- read.csv("goatiggc_cons_60k_11589869_A.csv", skip=7)

samples <- read.csv("Ziegen_HU-Berlin_samples.txt", sep="\t", row.names = 1)
samples[1:10,]

snpdata <- read.csv("Ziegen_HU-Berlin_Matrix.txt", sep="\t", skip = 9, header=TRUE, check.names=FALSE, row.names=1, na.strings=c("NA", "", "--"))
snpdata[1:10,1:10]

snpinfo <- read.csv("Ziegen_HU-Berlin_SNP.txt", sep="\t", row.names=1)
snpinfo[1:10,]
snpinfo <- snpinfo[which(snpinfo[,"GenTrain.Score"] >= 0.6),]             # GenTrain > 0.6 is a good SNP
snpinfo <- snpinfo[which(snpinfo[,"Minor.Freq"] >= 0.05),]                # Minimum allele frequency of minor allele 5 %
snpinfo <- snpinfo[-which(rownames(snpinfo) %in% posUncertain),]          # No uncertain positions
dim(snpinfo)

arrayinfo <- arrayinfo[match(rownames(snpinfo), arrayinfo[,2]),]
rownames(arrayinfo) <- arrayinfo[,2]
dim(arrayinfo)

ILMNannot <- ILMNannot[match(rownames(snpinfo), ILMNannot[,"Name"]),]
rownames(ILMNannot) <- ILMNannot[,"Name"]
dim(ILMNannot)

snpinfo[,"Chr"] <- arrayinfo[rownames(snpinfo),"chr"]
snpinfo[,"Position"] <- arrayinfo[rownames(snpinfo),"chr_pos"]
snpinfo <- cbind(snpinfo, allele = arrayinfo[rownames(snpinfo),"allele"])
snpinfo <- cbind(snpinfo, rsID = arrayinfo[rownames(snpinfo),"rs."])
snpinfo <- cbind(snpinfo, reference = arrayinfo[rownames(snpinfo),"reference"])
snpinfo <- cbind(snpinfo, SourceSeq = ILMNannot[rownames(snpinfo),"SourceSeq"])

snpinfo <- snpinfo[order(snpinfo[,"Chr"], snpinfo[,"Position"]),]
snpdata <- snpdata[rownames(snpinfo), ]


setwd("E:/Goat/DNA/Siham Analysis")
write.table(snpdata, "filtered_snps.txt", sep="\t", quote=FALSE)
write.table(snpinfo, "snpinfo.txt", sep="\t", quote=FALSE)
write.table(samples, "sampleinfo.txt", sep="\t", quote=FALSE)

