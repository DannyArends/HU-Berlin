source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")

# Load the genotype call data
phased.vcf <- read.table(gzfile(paste0("TRDredo/phased.vcf.gz")), header = FALSE, colClasses="character")     # Load
colnames(phased.vcf)  <- strsplit(sub("#","",readLines(gzfile(paste0("TRDredo/phased.vcf.gz")), n=10)[10]),"\t")[[1]]       # Add column header

samples <- colnames(phased.vcf)[-c(1:9)]

phased.vcf[, samples] <- apply(phased.vcf[, samples],2,function(x){unlist(lapply(strsplit(x,":"),"[",1))})                   # Only the genotypes
rownames(phased.vcf) <- phased.vcf[,"ID"]
phased.vcf <- phased.vcf[,-which(colnames(phased.vcf) %in% c("ID", "QUAL","FILTER","INFO","FORMAT"))]                        # Remove unneeded columns
phased.vcf[1:10,1:10]

# Load in the phenotypes and create the pedigree file
phenotypes <- read.table("Phenotypes/allPhenotypes.txt",sep="\t",header=TRUE,na.strings=c(0,"-","NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.vcf)),]                            # We do not have genotypes for all individuals

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals

# Change the genotype coding
phased.AHB <- phased.vcf                                                                         # Copy
phased.AHB[, samples] <- apply(phased.AHB[, samples], 2, fromVCFcode.Num)                        # Change coding to A H B

# QTL Analysis of generation 28
# Filter the genotypes for markers that are seggregating
phased.AHB <- phased.AHB[-which(lapply(apply(phased.AHB[, F2], 1, table),length) == 1),]

genotypes <- NULL
for(x in c(seq(1,19),"X","Y","M")) {  # Reorder the chromosomes and name it genotypes
  genotypes <- rbind(genotypes, phased.AHB[which(phased.AHB[,"CHROM"] == x), ])
}
numgeno <- genotypes[, -c(1:4)]
#numgeno <- read.table(file="genotypes_F2_filtered_ordered_numeric.txt", sep = '\t', check.names = FALSE)
numgeno[1:10,1:10]

Meff_PCA <- function(eigenValues, PCA_cutoff = 0.995){
  totalEigenValues <- sum(eigenValues)
  myCut <- PCA_cutoff*totalEigenValues
  num_Eigens <- length(eigenValues)
  myEigenSum <- 0
  index_Eigen <- 0
  
  for(i in 1:num_Eigens){
    if (myEigenSum <= myCut) {
      myEigenSum <- myEigenSum + eigenValues[i]
      index_Eigen <- i
    } else {
      break
    }
  }  
  return(index_Eigen)
}

inferCutoff <- function(dt_My, PCA_cutoff = 0.995){ # infer the cutoff => Meff
  CLD <- cor(dt_My, use = "pair")
  eigen_My <- eigen(CLD)
  eigenValues_dt <- abs(eigen_My$values)
  Meff_PCA_gao <- Meff_PCA(eigenValues_dt, PCA_cutoff)
  return(Meff_PCA_gao)
}

numLoci <- nrow(numgeno)
simpleMeffs <- NULL
rangeTested <- seq(100, 1200, 100)
for(x in rangeTested){
  simpleMeff <- NULL
  fixLength <- x 
  i <- 1
  myStart <- 1
  myStop <- 1
  while (myStop < numLoci) {
    myDiff <- numLoci - myStop 
    if (myDiff <= fixLength) break
    
    myStop <- myStart + i*fixLength - 1
    snpInBlk <- t(numgeno[myStart:myStop, ])
    MeffBlk <- inferCutoff(snpInBlk)
    simpleMeff <- c(simpleMeff, MeffBlk)
    myStart <- myStop+1
  }
  snpInBlk <- t(numgeno[myStart:numLoci, ])
  if (nrow(snpInBlk) > 1) {
    MeffBlk <- inferCutoff(snpInBlk)
  } else {
    MeffBlk <- 1
  }
  simpleMeff <- c(simpleMeff, MeffBlk)

  cat("Total number of SNPs is: ", numLoci, "\n")
  cat("Inferred Meff is: ", sum(simpleMeff), "\n")
  simpleMeffs <- c(simpleMeffs, sum(simpleMeff))
}
c(fixLength = rangeTested[which.min(simpleMeffs)], nTests = min(simpleMeffs))
plot(rangeTested, simpleMeffs, xlab = "fixLength", ylab="Meff")