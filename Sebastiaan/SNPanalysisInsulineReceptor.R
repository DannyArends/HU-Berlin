# SNPanalysisInsulineReceptor.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written May, 2014
#
# Analysis of data from Sebastiaan

setwd("E:/Atlas/")

nonDiabetic   <- c("BFMI852", "BFMI856", "BFMI860-12", "BFMI860-S2")                          # nonDiabetic BFMI mice
isDiabetic    <- c("BFMI861-S1")                                                              # Diabetic BFMI strain
isDuplicated  <- c("BFMI860-12", "BFMI861-S2", "BFMI861-S1")                                  # Strains measured multiple times

SNPdata <- read.table("SNPAnnotated.txt", sep="\t", header=TRUE)
cat("Starting with", nrow(SNPdata), "valid SNPs from the mouse diversity array\n")

annotation <- read.table("MouseAnnotation.txt", header=TRUE)                                  # Load the annotation
validAnnot <- match(colnames(SNPdata)[8:ncol(SNPdata)], colnames(annotation))                 # Annotation that matches the mice we have in our SNP data
annotation <- t(annotation[1:6, validAnnot])
colnames(annotation) <- c("Date", "Line", "Generation", "Phenotype", "Sex", "Diet")           # Correct header names

BFMIlines  <- grep("BFMI", annotation[,"Line"])                                               # Get only the BFMI lines
annotation <- annotation[BFMIlines, ]                                                         # Get only the BFMI lines

inconsistentSNPs <- NULL
for(mouseLine in isDuplicated){                                                               # Find the inconsistent SNPs between the strains
  mouseNames <- names(which(annotation[,"Line"] == mouseLine))
  unEqual <- which(apply(SNPdata[, mouseNames], 1, function(x){x[1] != x[2]}))
  cat("Found", length(unEqual), "differences between duplicates of", mouseLine,"\n")
  inconsistentSNPs <- c(inconsistentSNPs, unEqual)
}
inconsistentSNPs <- unique(inconsistentSNPs)

percentage <- paste0("(", round(length(inconsistentSNPs)/nrow(SNPdata)*100, d=1), "%)")       # % of inconsistent SNPs
cat("Removing", length(inconsistentSNPs), percentage, "inconsistently genotypes SNPs\n")
SNPdata <- SNPdata[-inconsistentSNPs, ]                                                       # Throw away the inconsistent SNPs

nonDMice <- NULL
for(mouseLine in nonDiabetic){
  nonDMice <- c(nonDMice, names(which(annotation[,"Line"] == mouseLine))[1])                  # The mouseID of a non-diabetic mouse
}

isDMouse <- names(which(annotation[,"Line"] == isDiabetic))[1]                                # The mouseID of the diabetic mouse

dSNP <- NULL                                                                                  # SNPs possibly involved in diabetic
for(snp in 1:nrow(SNPdata)){
  if(!(SNPdata[snp, isDMouse] %in% SNPdata[snp, nonDMice])){ dSNP <- c(dSNP, snp) }           # If the diabetic SNP is not in the nonDiabetic group, add the SNP
}
cat("Found", length(dSNP), "SNPs possibly involved with diabetes\n")

snpOUT <- cbind(SNPdata[dSNP,1:7], SNPdata[dSNP, isDMouse], SNPdata[dSNP, nonDMice])          # Create the ouput (subset the whole SNP arrays)
colnames(snpOUT) <- c(colnames(SNPdata)[1:7], isDiabetic, nonDiabetic)                        # Add a comprehensive header

NASNPinDiabetic <- which(is.na(snpOUT[,isDiabetic]))                                          # Not genotyped in the diabetic mouse                                        
cat("Removing", length(NASNPinDiabetic), "SNPs no genotype in the diabetic mouse\n")
snpOUT <- snpOUT[-NASNPinDiabetic,]                                                           # Remove them

NASNPinNonDiabetic <- which(apply(snpOUT[,nonDiabetic], 1, function(x){                       # Not genotyped in the non-diabetic mice
    sum(is.na(x))== length(nonDiabetic)})
  )
cat("Removing", length(NASNPinNonDiabetic), "SNPs no genotype in the non-diabetic mouse\n")
snpOUT <- snpOUT[-NASNPinNonDiabetic,]                                                        # Remoce them

percentage <- paste0("(", round(nrow(snpOUT) / nrow(SNPdata) * 100, d = 1), "%)")             # % of possibly diabetes causing SNPs
cat("Left with", nrow(snpOUT), percentage, "SNPs possibly involved in diabetes\n")
write.table(snpOUT, file="possibleDiabetesSNPs.txt", sep="\t", row.names = FALSE)             # Write them to a file

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

for(tid in 1:nrow(interestingTargets)){
  target   <- interestingTargets[tid,]
  chrSNPs  <- which(snpOUT[,"Chr"] == as.numeric(target[2]))
  snpOnCHR <- snpOUT[chrSNPs, ]
  inRange  <- ( snpOnCHR[, "Location"] > as.numeric(target[3]) - 100000 & snpOnCHR[, "Location"] < as.numeric(target[4]) + 100000 )
  if(any(inRange, na.rm=TRUE)){
    snpsInGene <- snpOnCHR[which(inRange), ]
    cat("Found a SNP near/in:",target[1],", SNP:", as.character(snpsInGene[,"dbSNP_ID"]),"\n")
  }
}

# Inside the gene regions +/- 100_000 basepairs:
#
# Found a SNP in: Cux1 , SNP: rs33588565  5:   136636006
# Found a SNP in: Pomc , SNP: rs46819789 12:     3985182 -> intron variant in EFR3
