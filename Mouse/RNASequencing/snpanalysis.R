# snpanalysis.R - Analyze the SNPs and indels called by the GenomeAnalysisToolKit
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

createNames <- function(x){ paste0(x[,1],":", x[,2],"_", x[,5]) }

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")

B6Nm <- read.table("Analysis/5068_GAGTGG_L004_.snps.vcf", colClasses="character")
B6Nf <- read.table("Analysis/5069_AGTCAA_L004_.snps.vcf", colClasses="character")

namesB6Nf <- createNames(B6Nf) ; namesB6Nm <- createNames(B6Nm)

B6N <- B6Nm[which(namesB6Nm %in% namesB6Nf),]
colnames(B6N) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

BFMIm <- read.table("Analysis/4868_GCCAAT_L001_.snps.vcf", colClasses="character")
BFMIf <- read.table("Analysis/5067_ATCACG_L004_.snps.vcf", colClasses="character")

namesBFMIf <- createNames(BFMIf) ; namesBFMIm <- createNames(BFMIm)

BFMI <- BFMIm[which(namesBFMIm %in% namesBFMIf),]
colnames(BFMI) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

# maternal B6N

matB6_1 <- read.table("Analysis/5070_CGATGT_L005_.snps.vcf", colClasses="character")
matB6_2 <- read.table("Analysis/5071_CCGTCC_L005_.snps.vcf", colClasses="character")
matB6_3 <- read.table("Analysis/5072_TAGCTT_L005_.snps.vcf", colClasses="character")

namesmB6_1 <- createNames(matB6_1) ; namesmB6_2 <- createNames(matB6_2) ; namesmB6_3 <- createNames(matB6_3)
matB6 <- matB6_1[which(namesmB6_1 %in% namesmB6_2 & namesmB6_1 %in% namesmB6_3),]
colnames(matB6) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

# maternal BFMI

matBFMI_1 <- read.table("Analysis/5073_TTAGGC_L006_.snps.vcf", colClasses="character")
matBFMI_2 <- read.table("Analysis/5074_GATCAG_L006_.snps.vcf", colClasses="character")
matBFMI_3 <- read.table("Analysis/5075_ATGTCA_L006_.snps.vcf", colClasses="character")

namesmBFMI_1 <- createNames(matBFMI_1) ; namesmBFMI_2 <- createNames(matBFMI_2) ; namesmBFMI_3 <- createNames(matBFMI_3)
matBFMI <- matB6_1[which(namesmBFMI_1 %in% namesmBFMI_2 & namesmBFMI_1 %in% namesmBFMI_3),]
colnames(matBFMI) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

rownames(B6N)   <- createNames(B6N)
rownames(BFMI)   <- createNames(BFMI)
rownames(matB6)   <- createNames(matB6)
rownames(matBFMI) <- createNames(matBFMI)

for(snp in rownames(matB6)){
  if(matB6[snp,"FORMAT"] == "GT:AD:DP:GQ:PL"){
    values <- strsplit(matB6[snp,"SAMPLE"], ":")
    totalReads <- as.numeric(unlist(values[[1]][3]))
    alleleReads <- as.numeric(unlist(strsplit(unlist(values[[1]][2]),",")))
    inBFMI  <- which(rownames(BFMI) ==  snp)
    inB6N   <- which(rownames(B6N) ==  snp)
    cat(snp,": ", totalReads,"->", alleleReads/totalReads, "\n")
  }
}

formatHighQual <- which(unique(matBFMI[,"FORMAT"]) == "GT:AD:DP:GQ:PL" & MIX[, "FILTER"] != "LowQual")

mmatrix <- NULL
for(x in 1:nrow(BFMI)){
  totalReads <-
  alleleReads <- as.numeric(unlist(strsplit(unlist(strsplit(BFMI[x,"SAMPLE"], ":")[[1]][2]),",")))
  cat(x,"has",totalReads,":", length(alleleReads), "=>", alleleReads/totalReads, "\n")
  mmatrix <- rbind(mmatrix, alleleReads/totalReads)
 # if(length(alleleReads) > 2) stop()
}


lapply(lapply(alleleReads, strsplit,","),function(x){ as.numeric(x[1]) / as.numeric(x[2])})

# Filer: B6N == B6N && BFMI == BFMI (not on X, Y)
# Filter first: B6N != BFMI
# Then: Filter 3-alleles (perhaps keep them in as 2 alleles when 1 is very low < 2 %)
# Then: Find the Frequencies of the remaining SNPs in F1 (by group)

