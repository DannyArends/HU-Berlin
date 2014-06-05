# SNPanalysisInsulineReceptor.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Analysis of data from Sebastiaan
# Insuline receptor: Chromosome 8: 3,150,922 - 3,279,617 reverse strand.

setwd("D://Sebastiaan_mouse")

chromosomes <- c(1:19, "M", "X", "Y")
rowheader   <- 1:6

lapply(chromosomes, function(x){
  chrdata <- read.csv(paste0("Mm_chr", x, "_snp_genotypes.txt"), colClasses = "character", sep="\t")

  chr <- chrdata[10:nrow(chrdata),]
  colnames(chr) <- as.character(unlist(chrdata[9,]))
  colnames(chr)[7:ncol(chr)] <- as.character(unlist(chrdata[4,7:ncol(chr)]))

  interesting <- c("BFMI861-S1","BFMI861-S2")
  snpchr <- chr[, interesting]

  dSNP   <- which(snpchr[,1]!=snpchr[,2])
  snpOut <- cbind(chr[dSNP, rowheader], snpchr[dSNP,])

  write.table(snpOut, file=paste0("MM_diffences_chr", x, "_snps.txt"), sep="\t",row.names = FALSE)
})

chr8 <- read.table("MM_diffences_chr8_snps.txt", header=TRUE)

snps <- which(chr8[,"SNP.Position..Build.37."] > 3150922 & chr8[,"SNP.Position..Build.37."] < 3279617)
snps   # We're in luck there is 1 SNP inside the insuling receptor, however its in an intron


snps <- which(chr8[,"SNP.Position..Build.37."] > 3279617 & chr8[,"SNP.Position..Build.37."] < 3303617)
snps   # We're in luck there is 1 SNP inside the insuling receptor, however its in an intron

# JAX00659248 C T rs50264941 8 3170638 CC 0


