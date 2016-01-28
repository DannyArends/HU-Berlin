#
# Analyse ATB in megaMuga data from multiple generations after beagle phasing
#

source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

# Load the genotype call data
setwd("E:/Mouse/DNA/MegaMuga/")
phased.vcf <- read.table(gzfile(paste0("Analysis/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("Analysis/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header

samples <- colnames(phased.vcf)[-c(1:9)]

phased.vcf[, samples] <- apply(phased.vcf[, samples],2,function(x){unlist(lapply(strsplit(x,":"),"[",1))})                   # Only the genotypes
rownames(phased.vcf) <- phased.vcf[,"ID"]
phased.vcf <- phased.vcf[,-which(colnames(phased.vcf) %in% c("ID", "QUAL","FILTER","INFO","FORMAT"))]                              # Remove unneeded columns
phased.vcf[1:10,1:10]

# Change the genotype coding
phased.AHBp <- phased.vcf                                                                         # Copy
phased.AHBp[, samples] <- apply(phased.AHBp[, samples], 2, fromVCFcode.AHBp)                      # Change coding to A H0 H1 B

# Change the genotype coding
phased.geno <- phased.vcf                                                                         # Copy
phased.geno[, samples] <- apply(phased.geno[, samples], 2, fromVCFcode.geno)                      # Change coding to AA AB BA BB


# Load the Allele Transmission Biased regions
sPat28 <- read.table("Analysis/TransmissionBias_Pat_0.01_28.txt", sep="\t")
sMat28 <- read.table("Analysis/TransmissionBias_Mat_0.01_28.txt", sep="\t")
# Load Phenotypes
phenotypes <- read.table("Phenotypes/allPhenotypes.txt", sep="\t", header=TRUE, na.strings=c(0, "-", "NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.AHBp)),]                # We do not have genotypes for all individuals
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                          # Add the season column to the matrix

individuals <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]

results <- NULL
for(m in rownames(sPat28)){
  g12 <- NULL
  g21 <- NULL
  for(i in individuals){
    vG <- phased.AHBp[m, as.character(phenotypes[i, "Vater"])]
    iG <- phased.AHBp[m, i]
    if(vG == "H0" || vG == "H1") {
      if(iG == "A" || iG == "H0") g21 <- c(g21, i)
      if(iG == "B" || iG == "H1") g12 <- c(g12, i)
    }
  }
  ind <- c(g12,g21)
  subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure  (factor)
  littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter   (linear effect)
  litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter     (factor)
  season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born     (factor)
  bfmiG          <- as.factor(unlist(phased.geno["UNC5048297", ind]))                                                       # Fixed effect: BFMI phased genotype (factor)
  genotype       <- as.factor(c(rep(1,length(g12)),rep(2, length(g21))))
  for(phe in c("d35", "d42", "d49", "d56", "d63", "d70","mri42d_fat","mri56d_fat","mri70d_fat")){
    res <- anova(lm(phenotypes[ind, phe] ~ littersize + litternumber + season + bfmiG + genotype)) 
    results <- rbind(results, c(m, phe, res[[5]]))
    cat(m, phe, res[[5]], "\n")
  }
  #cat(m, length(g12), length(g21),"\n")
}
colnames(results) <- c("Marker", "Phenotype", "LitterSize", "LitterNumber", "Season", "BFMIGenotype", "ATBgroup", "Residuals")
write.table(results, "Analysis/QTLmapping_PatRegions_BFMI_COF_NoSubFam.txt", sep ="\t", quote = FALSE, row.names = FALSE)

marker.annot <- read.table("Analysis/markerAnnotation.txt", sep="\t")
patQTL <- read.table("Analysis/QTLmapping_PatRegions_BFMI_COF_NoSubFam.txt", sep ="\t", row.names = NULL, header=TRUE)
matQTL <- read.table("Analysis/QTLmapping_MatRegions_BFMI_COF_NoSubFam.txt", sep ="\t", row.names = NULL, header=TRUE)

nPatRegions <- 49
nMatRegions <- 64

plotPheno <- function(patQTL, matQTL, phe = "d35"){
  op <- par(mfrow=c(2,1))
  iP <- which(patQTL[,"Phenotype"] == phe)
  plot(-log10(patQTL[iP,"ATBgroup"]), t ='h')
  iM <- which(matQTL[,"Phenotype"] == phe)
  plot(-log10(matQTL[iM,"ATBgroup"]), t ='h')
}

plotPheno(patQTL, matQTL)