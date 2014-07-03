# SNPanalysisDiabetes.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written May, 2014
#
# Analysis of data from Sebastiaan

setwd("E:/Mouse/DiversityArray/")

nonDiabetic   <- c("BFMI852", "BFMI856", "BFMI860-12", "BFMI860-S2")                              # nonDiabetic BFMI mice
isDiabetic    <- c("BFMI861-S1")                                                                  # Diabetic BFMI strain
isDuplicated  <- c("BFMI861-S2", "BFMI861-S1")                                                    # Strains measured multiple times

SNPdata <- read.table("Analysis/measurementsAtlas_annotated.txt", sep="\t", header=TRUE)
cat("Starting with", nrow(SNPdata), "valid SNPs from the mouse diversity array\n")

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                           # Load the annotation
validAnnot <- match(colnames(SNPdata)[9:ncol(SNPdata)], rownames(annotation))                     # Annotation that matches the mice we have in our SNP data
annotation <- annotation[validAnnot, ]

BFMIlines  <- grep("BFMI", annotation[,"Line"])                                                   # Get only the BFMI lines
annotation <- annotation[BFMIlines, ]                                                             # Get only the BFMI lines

inconsistentSNPs <- NULL
amount <- NULL
for(mouseLine in isDuplicated){                                                                   # Find the inconsistent SNPs between the strains
  mouseNames <- rownames(annotation[which(annotation[,"Line"] == mouseLine),])
  unEqual <- which(apply(SNPdata[, mouseNames], 1, function(x){x[1] != x[2]}))
  cat("Found", length(unEqual), "differences between duplicates of", mouseLine,"\n")
  amount <- c(amount, length(unEqual))
  inconsistentSNPs <- c(inconsistentSNPs, unEqual)
}
inconsistentSNPs <- unique(inconsistentSNPs)

percentage <- paste0("(", round(length(inconsistentSNPs)/nrow(SNPdata)*100, d=1), "%)")           # % of inconsistent SNPs
cat("Removing", length(inconsistentSNPs), percentage, "inconsistently genotyped SNPs\n")
SNPinconsistent <- SNPdata[inconsistentSNPs, ]                                                    # Inconsistent SNPs
write.table(SNPinconsistent[,1:8], file="Analysis/inconsistentSNPsAtlas.txt", sep="\t", row.names = FALSE)

dendrogram <- as.dendrogram(hclust(dist(t(SNPdata[,9:ncol(SNPdata)][,rownames(annotation)]))))    # Dendrogram of all BFMI lines
png(file="Analysis/Figures/DendrogramBFMI-lines.png")
  plot(dendrogram, main="Dendrogram BFMI-lines", xlab="", center=TRUE)
dev.off()
SNPdata <- SNPdata[-inconsistentSNPs, ]                                                           # Throw away the inconsistent SNPs

###### Comparison 1: BFMI861-S1 versus BFMI852, BFMI856, BFMI860-12, BFMI860-S2 ######

nonDMice <- NULL
for(mouseLine in nonDiabetic){
  nonDMice <- c(nonDMice, rownames(annotation[which(annotation[,"Line"] == mouseLine),])[1])      # The mouseID of a non-diabetic mouse
}

isDMouse <- rownames(annotation[which(annotation[,"Line"] == isDiabetic),])[1]                    # The mouseID of the diabetic mouse

dSNP <- NULL                                                                                      # SNPs possibly involved in diabetic
for(snp in 1:nrow(SNPdata)){
  if(!(SNPdata[snp, isDMouse] %in% SNPdata[snp, nonDMice])){ dSNP <- c(dSNP, snp) }               # If the diabetic SNP is not in the nonDiabetic group, add the SNP
}
cat("Found", length(dSNP), "SNPs possibly involved with diabetes\n")

Comparison1OUT <- cbind(SNPdata[dSNP,1:8], SNPdata[dSNP, isDMouse], SNPdata[dSNP, nonDMice])      # Create the ouput (subset the whole SNP arrays)
colnames(Comparison1OUT) <- c(colnames(SNPdata)[1:8], isDiabetic, nonDiabetic)                    # Add a comprehensive header

NASNPinDiabetic <- which(is.na(Comparison1OUT[,isDiabetic]))                                      # Not genotyped in the diabetic mouse                                        
cat("Removing", length(NASNPinDiabetic), "SNPs no genotype in the diabetic mouse\n")
Comparison1OUT <- Comparison1OUT[-NASNPinDiabetic,]                                               # Remove them

NASNPinNonDiabetic <- which(apply(Comparison1OUT[,nonDiabetic], 1, function(x){                   # Not genotyped in the non-diabetic mice
  sum(is.na(x))== length(nonDiabetic)
}))

cat("Removing", length(NASNPinNonDiabetic), "SNPs no genotype in the non-diabetic mouse\n")
Comparison1OUT <- Comparison1OUT[-NASNPinNonDiabetic,]                                            # Remove them

percentage <- paste0("(", round(nrow(Comparison1OUT) / nrow(SNPdata) * 100, d = 1), "%)")         # % of possibly diabetes causing SNPs
cat("Left with", nrow(Comparison1OUT), percentage, "SNPs possibly involved in diabetes\n")
write.table(Comparison1OUT, file="Analysis/Diabetes/BFMI861-S1vsALL_SNPs.txt", sep="\t", row.names = FALSE)

###### Comparison 2: BFMI861-S1 versus BFMI861-S2 ######

mouseS1 <- rownames(annotation[which(annotation[,"Line"] == "BFMI861-S1"),])[1]                   # The mouseID of the diabetic 861-S1 mouse
mouseS2 <- rownames(annotation[which(annotation[,"Line"] == "BFMI861-S2"),])[1]                   # The mouseID of the non-diabetic S2 mouse

dSNP <- NULL                                                                                      # SNPs different between S1 and S2
for(snp in 1:nrow(SNPdata)){
  if(!is.na(SNPdata[snp, mouseS1]) && !is.na(SNPdata[snp, mouseS2])){
    if(SNPdata[snp, mouseS1] != SNPdata[snp, mouseS2]){ dSNP <- c(dSNP, snp) }                    # If the SNP is different, add the SNP
  }
}
cat("Found", length(dSNP), "SNPs different between S1 and S2\n")

Comparison2OUT <- cbind(SNPdata[dSNP,1:8], SNPdata[dSNP, mouseS1], SNPdata[dSNP, mouseS2])        # Create the ouput (subset the whole SNP arrays)
colnames(Comparison2OUT) <- c(colnames(SNPdata)[1:8], "BFMI861-S1", "BFMI861-S2")                 # Add a comprehensive header
write.table(Comparison2OUT, file="Analysis/Diabetes/BFMI861-S1vsBFMI861-S2_SNPs.txt", sep="\t", row.names = FALSE)

###### Comparison 3: BFMI861-S1 equal to BFMI860-12 versus BFMI861-S2 ######

mouseS1 <- rownames(annotation[which(annotation[,"Line"] == "BFMI861-S1"),])[1]                   # The mouseID of the diabetic 861-S1 mouse
mouse86012 <- rownames(annotation[which(annotation[,"Line"] == "BFMI860-12"),])                   # The mouseIDs of the sensitive 860-12 mouse

nonSensitive  <- c("BFMI852", "BFMI856", "BFMI860-S2")                                            # nonSensitive BFMI mice
nonSMice <- NULL
for(mouseLine in nonSensitive){
  nonSMice <- c(nonSMice, rownames(annotation[which(annotation[,"Line"] == mouseLine),])[1])      # The mouseID of a non-diabetic mouse
}

dSNP <- NULL                                                                                      # SNPs different between S1 and S2
for(snp in 1:nrow(SNPdata)){
  if(!is.na(SNPdata[snp, mouseS1]) && !any(is.na(SNPdata[snp, mouse86012]))){
    if(SNPdata[snp, mouseS1] == SNPdata[snp, mouse86012[1]] && SNPdata[snp, mouseS1] == SNPdata[snp, mouse86012[2]]){      # BFMI861-S1 equal to BFMI860-12
      if(!(SNPdata[snp, mouseS1] %in% SNPdata[snp, nonSMice])){ dSNP <- c(dSNP, snp) }                                     # If the SNP is not in the nonSensitive group, add the SNP
    }
  }
}
cat("Found", length(dSNP), "SNPs equal between 861-S1 and 860-12, but different in other BFMI\n")

Comparison3OUT <- cbind(SNPdata[dSNP,1:8], SNPdata[dSNP, mouseS1], SNPdata[dSNP, mouse86012[1]], SNPdata[dSNP, nonSMice])  # Create the ouput (subset the whole SNP arrays)
colnames(Comparison3OUT) <- c(colnames(SNPdata)[1:8], "BFMI861-S1", "BFMI860-12", nonSensitive)                            # Add a comprehensive header
write.table(Comparison3OUT, file="Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_SNPs.txt", sep="\t", row.names = FALSE)
