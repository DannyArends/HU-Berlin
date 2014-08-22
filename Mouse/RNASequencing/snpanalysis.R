# snpanalysis.R - Analyze the SNPs and indels called by the GenomeAnalysisToolKit
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")

BFMI <- read.table("Analysis/4868_GCCAAT_L001_.snps.vcf", colClasses="character")
B6N <- read.table("Analysis/5068_GAGTGG_L004_.snps.vcf", colClasses="character")
MIX <- read.table("Analysis/5071.snps.recalibrated.vcf", colClasses="character")

BFMI[which(BFMI[,3] == "rs13463404"),]
which(B6N[,3] == "rs13463404")

colnames(MIX) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SAMPLE")

lowQuality <- which(MIX[, "FILTER"] == "LowQual")
