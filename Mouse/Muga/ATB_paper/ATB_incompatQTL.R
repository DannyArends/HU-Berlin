####################################################################
# Author    : Danny Arends
# Date      : 13-march-2016
# File Name : createMRI.R
# Purpose   : create an MRI table from the format provided by our MRI reader
# Used Files: dateToSeason.R
#             vcfTools.R
#############################################
#       Load external Functions             #
#############################################

source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

setwd("E:/Mouse/DNA/MegaMuga/")

# Load the Allele Transmission Biased regions
regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)

allregions <- rbind(regionsMat, regionsPat)

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
phased.geno <- phased.vcf                                                                         # Copy
phased.geno[, samples] <- apply(phased.geno[, samples], 2, fromVCFcode.AHB)                       # Change coding to AA AB BA BB

# Load Phenotypes
phenotypes <- read.table("Phenotypes/allPhenotypes.txt", sep="\t", header=TRUE, na.strings=c(0, "-", "NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- phenotypes[-which(!rownames(phenotypes) %in% colnames(phased.geno)),]                # We do not have genotypes for all individuals
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                          # Add the season column to the matrix

# Load the regions and their ChiSq scores
interacting <- read.table("ChiSqScores_Incompatible.txt",sep="\t")

phenoN <- c("d21", "d28", "d35", "d42", "mri42d_fat", "mri42d_lean", "d49", "d56", "mri56d_fat", "mri56d_lean", "d63", "d70", "mri70d_fat", "mri70d_lean")

# Get the BFMI individuals
ind <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]
bfmiG <- as.factor(unlist(phased.geno["UNC5048297", ind]))                                                                # BFMI phased genotype
ind <- ind[which(bfmiG == "H")]

subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure  (factor)
littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter   (linear effect)
litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter     (factor)
season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born     (factor)

res1 <-  vector("list", length(phenoN))
names(res1) <- phenoN
cnt <- 1
for(phe in phenoN) {
  for(x in 1:(nrow(allregions)-1)) {
    nm1 <- as.character(allregions[x,1])
    m1 <- as.factor(unlist(phased.geno[nm1,ind]))
    for(y in (x+1):nrow(allregions)) {
      nm2 <- as.character(allregions[y,1])
      m2 <- as.factor(unlist(phased.geno[nm2,ind]))
      mymodel <- lm(phenotypes[ind, phe] ~ subfamily + litternumber + m1*m2)
      myanova <- anova(mymodel)
      myAIC <- AIC(mymodel)
      #cat(nm1, nm2, phe, "LOD: ", round(unlist(-log10(myanova[[5]])), 2),"\n")
      #cat(nm1, nm2, phe, "VAR: ", round(100* (myanova[[2]] / sum(myanova[[2]])), 1), "\n")
      #cat(nm1, nm2, phe, "AIC: ", myAIC, "\n")
      resline <- c(nm1,nm2, round(unlist(-log10(myanova[[5]])), 2), myAIC)
      if(length(resline) == 9) res1[[cnt]] <- rbind(res1[[cnt]], resline)
    }
  }
  cat("Finished", phe, "\n")
  cnt <- cnt + 1
}

#nm1 <- "UNC30537625"
#nm2 <- "UNC18371164"
#m1 <- as.factor(unlist(phased.geno[nm1,ind]))
#m2 <- as.factor(unlist(phased.geno[nm2,ind]))
#mymodel <- lm(phenotypes[ind, phe] ~ subfamily + littersize + season + m1*m2)
#myanova <- anova(mymodel)
#round(100* (myanova[[2]] / sum(myanova[[2]])), 1)
#myAIC <- AIC(mymodel)


# Get the non BFMI individuals
ind <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]
bfmiG <- as.factor(unlist(phased.geno["UNC5048297", ind]))                                                       # Fixed effect: BFMI phased genotype (factor)
ind <- ind[which(bfmiG != "A")]

subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure  (factor)
littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter   (linear effect)
litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter     (factor)
season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born     (factor)

res2 <-  vector("list", length(phenoN))
names(res2) <- phenoN
cnt <- 1
for(phe in phenoN){
  for(x in 1:(nrow(allregions)-1)) {
    nm1 <- as.character(allregions[x,1])
    m1 <- as.factor(unlist(phased.geno[nm1,ind]))
    for(y in (x+1):nrow(allregions)) {
      nm2 <- as.character(allregions[y,1])
      m2 <- as.factor(unlist(phased.geno[nm2,ind]))
      mymodel <- lm(phenotypes[ind, phe] ~ subfamily + litternumber + m1*m2)
      myanova <- anova(mymodel)
      myAIC <- AIC(mymodel)
      #cat(nm1, nm2, phe, "LOD: ", round(unlist(-log10(myanova[[5]])), 2),"\n")
      #cat(nm1, nm2, phe, "VAR: ", round(100* (myanova[[2]] / sum(myanova[[2]])), 1), "\n")
      #cat(nm1, nm2, phe, "AIC: ", myAIC, "\n")
      resline <- c(nm1,nm2, round(unlist(-log10(myanova[[5]])), 2), myAIC)
      if(length(resline) == 9) res2[[cnt]] <- rbind(res2[[cnt]], resline)
    }
  }
  cat("Finished", phe, "\n")
  cnt <- cnt + 1
}

save(res1, file="TRD_pairwiseQTL_BFMI.RData")
save(res2, file="TRD_pairwiseQTL_nonBFMI.RData")

sign1 <- lapply(res1, function(x){
  threshold <- -log10(0.5/nrow(x))
  ii <- which(as.numeric(x[,8]) >= threshold)
  return(x[ii,])
})

sign2 <- lapply(res2, function(x){
  threshold <- -log10(0.7/nrow(x))
  ii <- which(as.numeric(x[,8]) >= threshold)
  return(x[ii,])
})

resNonBFMI <- NULL
for(x in 1:length(sign2)){
  if(length(sign2[[x]]) > 10){
    for(y in 1:nrow(sign2[[x]])){
      resNonBFMI <- rbind(resNonBFMI, c(names(sign2)[x], sign2[[x]][y,]))
    }
  }
  if(length(sign2[[x]]) == 10){
      resNonBFMI <- rbind(resNonBFMI, c(names(sign2)[x], sign2[[x]]))
  }
}

resBFMI <- NULL
for(x in 1:nrow(resNonBFMI)){
  subsetX <- res1[[resNonBFMI[x,1]]]
  resBFMI <- rbind(resBFMI, c(resNonBFMI[x,1], subsetX[which(subsetX[,1] == resNonBFMI[x,2] & subsetX[,2] == resNonBFMI[x,3]),]))
}

phased.geno[unique(c(resNonBFMI[,2],resNonBFMI[,3])),1:4]

# TODO Check the models, see if they make sense

nm1 <-"JAX00168530"
nm2 <-"UNC14860494"
m1 <- as.factor(unlist(phased.geno[nm1,ind]))
m2 <- as.factor(unlist(phased.geno[nm2,ind]))
mymodel <- lm(phenotypes[ind, phe] ~ subfamily + littersize + season + m1*m2)
myanova <- anova(mymodel)
round(100* (myanova[[2]] / sum(myanova[[2]])), 1)
myAIC <- AIC(mymodel)























mo <- lm(phenotypes[ind, phe] ~ subfamily)
boxplot((mo$residuals + mo$coefficients["(Intercept)"]) ~ subfamily)


m <- round(runif(length(ind)),0)

mo <- lm(phenotypes[ind, phe] ~ m)
boxplot(phenotypes[ind, phe] ~ m)
boxplot((mo$residuals + mo$coefficients["(Intercept)"]) ~ m)



gen28 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]
bfmiG <- as.character(unlist(phased.geno["UNC5048297", ind]))

gen27f <- as.character(phenotypes[ind, "Vater"])
bfmiGF <- as.character(unlist(phased.geno["UNC5048297", gen27f]))

gen27m <- as.character(phenotypes[ind, "Mutter"])
bfmiGM <- as.character(unlist(phased.geno["UNC5048297", gen27m]))

subFam <- as.factor(paste0(bfmiGF, bfmiGM, bfmiG))

subFam <- as.factor(paste0(bfmiGF, bfmiGM))


mo <- lm(phenotypes[gen28, phe] ~ subFam + bfmiG + litternumber)
myanova <- anova(mo)
round(100* (myanova[[2]] / sum(myanova[[2]])), 1)
AIC(mo)


boxplot(phenotypes[gen28, phe] ~ subFam)
boxplot((mo$residuals + mo$coefficients["(Intercept)"]) ~ subFam)
















# Get the non BFMI individuals
ind <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]
bfmiG <- as.factor(unlist(phased.geno["UNC5048297", ind]))                                                       # Fixed effect: BFMI phased genotype (factor)

subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure  (factor)
littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter   (linear effect)
litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter     (factor)
season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born     (factor)

res3 <-  vector("list", length(phenoN))
names(res3) <- phenoN
cnt <- 1
for(phe in phenoN){
  for(x in 1:(nrow(allregions)-1)) {
    nm1 <- as.character(allregions[x,1])
    m1 <- as.factor(unlist(phased.geno[nm1,ind]))
    for(y in (x+1):nrow(allregions)) {
      nm2 <- as.character(allregions[y,1])
      m2 <- as.factor(unlist(phased.geno[nm2,ind]))
      mymodel <- lm(phenotypes[ind, phe] ~ subfamily + bfmiG + litternumber + m1*m2)
      myanova <- anova(mymodel)
      myAIC <- AIC(mymodel)
      #cat(nm1, nm2, phe, "LOD: ", round(unlist(-log10(myanova[[5]])), 2),"\n")
      #cat(nm1, nm2, phe, "VAR: ", round(100* (myanova[[2]] / sum(myanova[[2]])), 1), "\n")
      #cat(nm1, nm2, phe, "AIC: ", myAIC, "\n")
      resline <- c(nm1,nm2, round(unlist(-log10(myanova[[5]])), 2), myAIC)
      if(length(resline) == 10) res3[[cnt]] <- rbind(res3[[cnt]], resline)
    }
  }
  cat("Finished", phe, "\n")
  cnt <- cnt + 1
}





