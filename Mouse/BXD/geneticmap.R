library(BXDtools)
setwd("D:/Edrive/Mouse/BxD/Juli2018")
bxd.genotypes <- read.table("bxd_geno_18Jan.txt", sep = "\t",skip=6, header=TRUE, na.strings=c("NA", "", "U"), colClasses="character")
bxd.genotypes <- bxd.genotypes[,-1]

physicalMap <- bxd.genotypes[,c(1, 2, 4, 5)]
colnames(physicalMap)[4] <- "cM"

bxd.genotypes <- bxd.genotypes[,-c(1:15)]
verbose <- FALSE
count.Heterozygous <- FALSE
start.at.zero <- FALSE

attr(bxd.genotypes, "map") <- physicalMap

geneticMap <- calculate.cM.positions(bxd.genotypes, FALSE, FALSE, TRUE)
