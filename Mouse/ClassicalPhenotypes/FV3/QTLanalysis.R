# QTL analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written March, 2009
#

library(qtl)
setwd("E:/Mouse/ClassicalPhenotypes/FV3")

# crossFF <- read.cross("csv", file="cross_F2_FF.csv",genotypes=c("A","H","B"), na.strings="NA")

crossNF <- read.cross("csv", file="cross_F2_NF.csv",genotypes=c("A","H","B"), na.strings="NA")
crossNF <- calc.genoprob(crossNF)

sex <- as.numeric(crossNF$pheno[,"Sex"])
season <- as.numeric(crossNF$pheno[,"sea"])
littersize <- as.numeric(crossNF$pheno[,"WG21"])

crossNF$pheno <- cbind(crossNF$pheno, FATDLEAN = crossNF$pheno[,"FAT70"] / crossNF$pheno[,"LEAN70"])
crossNF$pheno <- cbind(crossNF$pheno, FATDCW = crossNF$pheno[,"FAT70"] / crossNF$pheno[,"BW70_10"])

resFAT        <- scanone(crossNF, pheno.col="FAT70", intcovar=cbind(sex, season, littersize))
resFATDLEAN   <- scanone(crossNF, pheno.col="FATDLEAN", intcovar=cbind(sex, season, littersize))
resFATDCW     <- scanone(crossNF, pheno.col="FATDCW", intcovar=cbind(sex, season, littersize))
plot(resFATDLEAN, resFAT, resFATDCW)

RESvar <- lm(crossNF$pheno[,"FAT70"] ~ pull.geno(crossNF)[,"rs3151604"])

crossNF$pheno <- cbind(crossNF$pheno, ResFATDCW = RESvar / crossNF$pheno[,"CW"])

#SSexp   <- anova(lm(crossNF$pheno[,"FAT70"] ~ as.factor(pull.geno(crossNF)[,"rs3151604"])))[,"Sum Sq"]
SSexp   <- anova(lm(crossNF$pheno[,"FAT70"] ~ pull.geno(crossNF)[,"rs3151604"]))[,"Sum Sq"]
SStotal <- sum((crossNF$pheno[,"FAT70"] - mean(crossNF$pheno[,"FAT70"],na.rm=TRUE))^2,na.rm=TRUE)
SSexp/SStotal