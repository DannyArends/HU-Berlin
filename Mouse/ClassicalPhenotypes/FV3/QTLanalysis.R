# QTL analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written March, 2009
#

library(qtl)
setwd("E:/Mouse/ClassicalPhenotypes/FV3")

cross <- read.cross("csv", file="cross_F2.csv",genotypes=c("A","H","B"), na.strings="NA")
cross <- jittermap(cross)

cross$pheno <- cbind(cross$pheno, FATDLEAN = cross$pheno[,"FAT70"] / cross$pheno[,"LEAN70"])

sex         <- as.numeric(cross$pheno[,"Sex"])
season      <- as.numeric(cross$pheno[,"sea"])
futter      <- as.numeric(cross$pheno[,"Futter"])

resFATDLEANFUTTER <- scanone(cross, pheno.col="FATDLEAN", addcovar = cbind(sex, season, futter), intcovar = futter)
resFATDLEAN       <- scanone(cross, pheno.col="FATDLEAN", addcovar = cbind(sex, season, futter))
resFUTTER         <- scanone(cross, pheno.col="FATDLEAN", addcovar = cbind(sex, season))

plot(resFUTTER, resFATDLEAN, resFATDLEANFUTTER, main="Fat/Lean = Sex + Season + Futter + G + G:Futter")



crossFF <- read.cross("csv", file="cross_F2_FF.csv",genotypes=c("A","H","B"), na.strings="NA")
crossNF <- read.cross("csv", file="cross_F2_NF.csv",genotypes=c("A","H","B"), na.strings="NA")
crossNF <- calc.genoprob(crossNF)

crossNF$pheno <- cbind(crossNF$pheno, FATDLEAN = crossNF$pheno[,"FAT70"] / crossNF$pheno[,"LEAN70"])

genotypes <- pull.geno(fill.geno(crossNF))
phenotype <- crossNF$pheno[,"FATDLEAN"]

sex <- as.numeric(crossNF$pheno[,"Sex"])
season <- as.numeric(crossNF$pheno[,"sea"])
littersize <- as.numeric(crossNF$pheno[,"WG21"])

resFATDLEAN   <- scanone(crossNF, pheno.col="FATDLEAN", addcovar=cbind(sex, season, littersize))
topmarker <- genotypes[,rownames(resFATDLEAN[which.max(resFATDLEAN[,3]),])]

genotypes <- genotypes[which(topmarker!=1),]
phenotype <- phenotype[which(topmarker!=1)]

lods <- apply(genotypes, 2, function(genotype){ return(-log10(anova(lm(phenotype ~ as.factor(genotype)))[[5]][1])) })
plot(lods, t='l')

crossNF$pheno <- cbind(crossNF$pheno, FATDCW = crossNF$pheno[,"FAT70"] / crossNF$pheno[,"CW"])

resFAT        <- scanone(crossNF, pheno.col="FAT70", addcovar=cbind(sex, season, littersize))

plot(resFAT, resFATDLEAN, col=c("green", "blue"))

RESvar <- lm(crossNF$pheno[,"FATDLEAN"] ~ sex + season + littersize + genotypes[,which(resFATDLEAN[,3] > 10)][,1])
phenoresiduals <- rep(NA,nrow(crossNF$pheno))
phenoresiduals[as.numeric(names(RESvar$residuals))] <- RESvar$residuals

crossNF$pheno   <- cbind(crossNF$pheno, FATDLEANres = phenoresiduals)
resFATDLEANres  <- scanone(crossNF, pheno.col="FATDLEANres")
plot(resFATDLEAN, resFATDLEANres, col=c("blue", "black"))


SSexp   <- anova(lm(crossNF$pheno[,"FATDLEANres"] ~ as.factor(pull.geno(crossNF)[,"rs3151604"])))[,"Sum Sq"]
SStotal <- sum((crossNF$pheno[,"FATDLEANres"] - mean(crossNF$pheno[,"FATDLEANres"],na.rm=TRUE))^2,na.rm=TRUE)
SSexp/SStotal


