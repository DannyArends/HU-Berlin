# Preprocessing of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016

setwd("D:/Edrive/Goat/DNA/Annotation")

arrayinfo <- read.csv("IGGC CapriBatch1.txt", sep="\t", skip = 29)
arrayinfo[1:10,]
posUncertain <- as.character(arrayinfo[which(duplicated(arrayinfo[,"loc_snp_id"])),2])

arrayinfo <- cbind(arrayinfo, reference = NA)

for(chr in c(1:29, "X")){
  dnasequence <- readLines(gzfile(paste0("CHIR1.0/chr",chr,".fna.gz")))       # Read the per chromosome fasta file
  dnasequence <- paste0(dnasequence[-1], collapse="")                         # Remove the first line, and collapse the rest of the data
  dnasequence <- strsplit(dnasequence, "")[[1]]                               # Split the DNA string into individual letters
  arrayinfo[arrayinfo[,"chr"] == chr, "reference"] <- dnasequence[arrayinfo[arrayinfo[,"chr"] == chr,"chr_pos"]]
}

setwd("E:/Goat/DNA/SihamRawData")

ILMNannot <- read.csv("goatiggc_cons_60k_11589869_A.csv", skip=7)

samples <- read.csv("Ziegen_HU-Berlin_samples.txt", sep="\t", row.names = 1)
samples[1:10,]

snpdata <- read.csv("Ziegen_HU-Berlin_Matrix.txt", sep="\t", skip = 9, header=TRUE, check.names=FALSE, row.names=1, na.strings=c("NA", "", "--"))
snpdata[1:10,1:10]

snpinfo <- read.csv("Ziegen_HU-Berlin_SNP.txt", sep="\t", row.names=1)
snpinfo[1:10,]

dim(snpinfo)
snpinfo <- snpinfo[which(snpinfo[,"GenTrain.Score"] >= 0.6),]             # GenTrain >= 0.6 is a good SNP
dim(snpinfo)
snpinfo <- snpinfo[which(snpinfo[,"Minor.Freq"] >= 0.05),]                # Minimum allele frequency of minor allele 5 %
dim(snpinfo)
snpinfo <- snpinfo[-which(rownames(snpinfo) %in% posUncertain),]          # No multi mapping positions
dim(snpinfo)

# Create the annotation files, based on the current number of good SNPs
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

# We still need to get rid of SNPs with > 5% missing data
callrates <- (apply(apply(snpdata,1,is.na),2,sum) / ncol(snpdata))                # Call rates
snpdata <- snpdata[-which(callrates > 0.05),]                                     # Require a maximum of 5 % missing data at a marker
snpinfo <- snpinfo[rownames(snpdata),]
dim(snpdata)

noPos <- which(is.na(snpinfo[,"Position"]))                                       # unknown posiiton (remove them as well)
snpinfo <- snpinfo[-noPos,]
snpdata <- snpdata[rownames(snpinfo), ]
dim(snpdata)


setwd("E:/Goat/DNA/SihamAnalysis")
write.table(snpdata, "filtered_snps.txt", sep="\t", quote=FALSE)
write.table(snpinfo, "snpinfo.txt", sep="\t", quote=FALSE)
write.table(samples, "sampleinfo.txt", sep="\t", quote=FALSE)




setwd("E:/Goat/DNA/SihamAnalysis")
snpdata <- read.table("filtered_snps.txt", sep="\t", check.names=FALSE)
snpdata <- snpdata[, -which(colnames(snpdata) == "DN 2")]
write.table(snpdata, "filtered_snps_NO_DN2.txt", sep="\t", quote=FALSE)