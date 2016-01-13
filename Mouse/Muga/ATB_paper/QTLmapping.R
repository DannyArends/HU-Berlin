
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

# Load in the phenotypes and create the pedigree file
phenotypes <- read.table("Phenotypes/allPhenotypes.txt",sep="\t",header=TRUE,na.strings=c(0,"-","NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.vcf)),]                # We do not have genotypes for all individuals

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals

# Change the genotype coding
phased.AHB <- phased.vcf                                                                         # Copy
phased.AHB[, samples] <- apply(phased.AHB[, samples], 2, fromVCFcode.AHB)                        # Change coding to A H B


# QTL Analysis of generation 28
# Filter the genotypes for markers that are seggregating
phased.AHB <- phased.AHB[-which(lapply(apply(phased.AHB[, F2], 1, table),length) == 1),]

genotypes <- NULL
for(x in c(seq(1,19),"X","Y","M")) {  # Reorder the chromosomes and name it genotypes
  genotypes <- rbind(genotypes, phased.AHB[which(phased.AHB[,"CHROM"] == x), ])
}

# Now we want to do the QTL mapping for a single phenotype
covariates <- phenotypes[F2, c("Eltern_ID", "WG2", "W.Label", "Season")]
phenotype  <- phenotypes[F2, "mri70d_fat"]

getParent <- function(genotypes, phenotypes, individuals = F2, marker = "UNC5048297", parent = "Mutter"){
  genotypes[marker, as.character(phenotypes[F2, parent])]
}

parent <- getParent(genotypes, phenotypes)

cnt <- 1
pvalues <- t(apply(genotypes[, F2], 1, function(marker){
  mylm <- lm(phenotype ~ covariates[,"Eltern_ID"] + covariates[,"WG2"] + covariates[,"W.Label"] + covariates[,"Season"] + as.factor(unlist(parent)) * as.factor(as.character(marker)))
  cnt <<- cnt + 1
  r <- unlist(anova(mylm)[[5]])
  if(length(r) == 7) return(r)
  return(rep(NA,7))
}))

plot(-log10(p.adjust(pvalues[,6],"BH")), col = as.numeric(as.factor(genotypes[,"CHROM"])),t ='h')
points(-log10(p.adjust(pvalues[,5],"BH")), col = as.numeric(as.factor(genotypes[,"CHROM"])),t ='l')
points(-log10(p.adjust(pvalues[,7],"BH")), col = as.numeric(as.factor(genotypes[,"CHROM"])),t ='l')
