# Test the 5 different POE models at the chromosome 3 locus
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sept, 2015
# first written Sept, 2015


source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypesFP <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("", "??"))                            # Phased by Beagle (ALL)
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes

missingPerInd <- apply(genotypesFP, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("33310233", "6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                         # Remove the faulty individuals from the phenotypes
genotypesFP <- genotypesFP[,-which(colnames(genotypesFP) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes


# We can test for transmission bias: e.g. is the father more likely to transmit the BFMI allele ?
for(x in 1:100){
  ch <- genotypesFP[x,"6661223"]
  vat <- genotypes[x,as.character(phenotypes["6661223","Vater"])]
  mut <- genotypes[x,as.character(phenotypes["6661223","Mutter"])]
  cat(ch,vat,mut,"\n")
}



F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]; length(F2)                                     # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]; length(F1)                                     # The F1 individuals
P <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]; length(P)                                       # The P individuals

mothers <- F1[which(phenotypes[F1,"sex"] == "f")]
fathers <- F1[which(phenotypes[F1,"sex"] == "m")]

topM <- "UNC5048297"


mA <- colnames(genotypes[,mothers])[which(genotypes[topM,mothers] == "A")]
mB <- colnames(genotypes[,mothers])[which(genotypes[topM,mothers] == "B")]

fA <- colnames(genotypes[,fathers])[which(genotypes[topM,fathers] == "A")]
fB <- colnames(genotypes[,fathers])[which(genotypes[topM,fathers] == "B")]

cfA <- F2[which(phenotypes[F2,c("Vater")] %in% fA)]
cfB <- F2[which(phenotypes[F2,c("Vater")] %in% fB)]

cmA <- F2[which(phenotypes[F2,c("Mutter")] %in% mA)]
cmB <- F2[which(phenotypes[F2,c("Mutter")] %in% mB)]

cmB %in% cfA


