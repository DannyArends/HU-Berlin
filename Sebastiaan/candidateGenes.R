# candidateGenes.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of candidate gene data from Sebastiaan

# Might be interesting as targets, GRCm38.p2 coordinates
#
# Insr: Insuline receptor       Chromosome 8:    3,150,922 -   3,279,617 reverse strand.
# GLUT4/Slc2a4 :                Chromosome 11:  69,942,539 -  69,948,188
#
# From the article: Stratigopoulos et al. (Cell metabolism)
# 
# FTO:                          Chromosome 8:   91,313,532 -  91,668,439
# Rpgrip1l:                     Chromosome 8:   91,217,030 -  91,313,262 reverse strand
# Cux1:                         Chromosome 5:  136,248,135 - 136,567,490 reverse strand
# Lepr:                         Chromosome 4:  101,717,404 - 101,815,352
# STAT3:                        Chromosome 11: 100,885,098 - 100,939,540 reverse strand
#
# Npy:                          Chromosome 6:   49,822,710 -  49,829,507
# Agrp:                         Chromosome 8:  105,566,700 - 105,568,298 reverse strand
# Pomc:                         Chromosome 12:   3,954,951 -   3,960,618
# Mc4r:                         Chromosome 18:  66,857,715 -  66,860,472 reverse strand

interestingTargets <- rbind(c("Insr",          8,   3150922,   3279617),
                            c("GLUT4/Slc2a4", 11,  69942539,  69948188),
                            c("FTO",           8,  91313532,  91668439),
                            c("Rpgrip1l",      8,  91217030,  91313262),
                            c("Cux1",          5, 136248135, 136567490),
                            c("Lepr",          4, 101717404, 101815352),
                            c("STAT3",        11, 100885098, 100939540),
                            c("Npy",           6,  49822710,  49829507),
                            c("Agrp",          8, 105566700, 105568298),
                            c("Pomc",         12,   3954951,   3960618),
                            c("Mc4r",         18,  66857715,  66860472))

setwd("E:/Atlas/")
snpOUT       <- read.table(file, sep="\t", header=TRUE)

for(tid in 1:nrow(interestingTargets)){
  target   <- interestingTargets[tid,]
  chrSNPs  <- which(snpOUT[,"Chr"] == as.numeric(target[2]))
  snpOnCHR <- snpOUT[chrSNPs, ]
  inRange  <- ( snpOnCHR[, "Location"] > as.numeric(target[3]) - 20000 & snpOnCHR[, "Location"] < as.numeric(target[4]) + 20000 )
  if(any(inRange, na.rm=TRUE)){
    snpsInGene <- snpOnCHR[which(inRange), ]
    cat("Found a SNP near/in:",target[1],", SNP:", as.character(snpsInGene[,"dbSNP_ID"]),"\n")
  }
}
