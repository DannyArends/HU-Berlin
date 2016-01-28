
source("D:/Github/HU-Berlin/Mouse/Muga/ATB_paper/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes
genotypes2   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes
genotypesPh <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("", "AA", "CC", "TT", "GG", "??"))    # Phased by Beagle (only heterozygous)
genotypesFP <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("", "??"))                            # Phased by Beagle (ALL)
genotypesGP <- read.table("Analysis/genotypesPhasedGP.txt",     sep="\t", check.names=FALSE, colClasses="character", na.strings=c(""))                                  # Phased towards the grandparents

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("33310233", "6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                         # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]; length(F2)                                     # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]; length(F1)                                     # The F1 individuals
P <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]; length(P)                                       # The P individuals

## Lecture
#subsetMkr <- sample(rownames(subsetG),3000)

#subsetP <- phenotypes[F2, c(paste0("d",seq(21,70,7)),"WG2")]
#subsetG <- genotypes[subsetMkr, F2]
#subsetM <- map[subsetMkr,]
#setwd("D:/Projects/Lectures 2014-2015/R course/R course/Assignments")

#write.table(subsetP, "lecture5_phenotypes.txt", sep="\t", quote=FALSE)
#write.table(subsetG, "lecture5_genotypes.txt", sep="\t", quote=FALSE)
#write.table(subsetM, "lecture5_map.txt", sep="\t", quote=FALSE)

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
birthmonth <- unlist(lapply(strsplit(as.character(phenotypes[,"W.dat"]),".", fixed=TRUE),"[",2))
phenotypes <- cbind(phenotypes, Birthmonth = birthmonth)                                                      # Add the birth month column to the matrix

phenotypes <- cbind(phenotypes, mri42d_fatDlean = phenotypes[,"mri42d_fat"] / phenotypes[,"mri42d_lean"])     # Fat / Lean day 42
phenotypes <- cbind(phenotypes, mri56d_fatDlean = phenotypes[,"mri56d_fat"] / phenotypes[,"mri56d_lean"])     # Fat / Lean day 56
phenotypes <- cbind(phenotypes, mri70d_fatDlean = phenotypes[,"mri70d_fat"] / phenotypes[,"mri70d_lean"])     # Fat / Lean day 70

onegenotype <- which(lapply(apply(genotypes[,F2], 1, table), length) < 2)                                       # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- (genotypes[-onegenotype, F2])                                                                    # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypes), "markers\n")                                                                  # == Left with 21260 markers    // Left with 11581 markers

onegenotype <- which(lapply(apply(genotypesPh[,F2], 1, table), length) < 2)                                     # Markers with only one genotype cannot be used in QTL mapping
genotypesPh <- (genotypesPh[-onegenotype,F2])                                                                   # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesPh), "markers\n")                                                                # == Left with 20199 markers    // Left with 19465 markers

onegenotype <- which(lapply(apply(genotypesFP[,F2], 1, table), length) < 2)                                     # Markers with only one genotype cannot be used in QTL mapping
genotypesFP <- (genotypesFP[-onegenotype,F2])                                                                   # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesFP), "markers\n")                                                                # == Left with 20930 markers    // Left with 19908 markers

onegenotype <- which(lapply(apply(genotypesGP[,F2], 1, table), length) < 2)                                     # Markers with only one genotype cannot be used in QTL mapping
genotypesGP <- (genotypesGP[-onegenotype,F2])                                                                   # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypesGP), "markers\n")                                                                # == Left with 19903 markers    // Left with 7551 markers

phenotypes <- phenotypes[F2,]

#genotypes <- genotypes[, which(genotypes["UNC4784019", F2] != "A")]
onX <- rownames(map[which(map[,1] == "X"),])
onX <- onX[which(onX %in% rownames(genotypes))]

results <- NULL
for(x in onX){
  ind            <- colnames(genotypes[x,!is.na(genotypes[x,])])                                                            # Which individuals have genotype data
  subfamily      <- as.factor(phenotypes[ind, "Vater"])                                                                     # Fixed effect: Subfamily structure (factor)
  littersize     <- as.numeric(phenotypes[ind, "WG2"])                                                                      # Fixed effect: Size of the litter  (linear effect)
  litternumber   <- as.factor(phenotypes[ind, "W.Label"])                                                                   # Fixed effect: Number of litter    (factor)
  season         <- as.factor(phenotypes[ind, "Season"])                                                                    # Fixed effect: Season when born    (factor)
  genotype       <- as.factor(t(genotypes[x,!is.na(genotypes[x,])]))                                                        # The genotype under investigation  (factor)

  phenotype      <- phenotypes[ind, "d70"]                      #/ phenotypes[ind, paste0("mri",pheno.col,"_lean")]     # Response: Fat / Lean
  tryCatch(
        res <- anova(lm(phenotype ~ littersize + litternumber + genotype))[[5]], 
    error = function(e){ res <<- rep(NA, 2) })
  results <- rbind(results, res)
}

rownames(results) <- onX


 cbind(map[names(results[which(-log10(results[,3]) > 2.5),3]),],-log10(results[which(-log10(results[,3]) > 2.5),3]))