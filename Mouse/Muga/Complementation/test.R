####################################################################
# Author    : Danny Arends
# Date      : 13-march-2016
# File Name : test.R
# Purpose   : 
# Used Files: createMRI.R
#############################################
#              packages & external          #
#############################################
source("D:\\Ddrive\\Github\\HU-Berlin\\Mouse\\Muga\\Complementation\\createMRI.R")

#############################################
#                Functions                  #
#############################################

adjust.phenotype <- function(geno, tp = "21"){
  model <-  lm(geno[,tp] ~ as.factor(geno[,"WG"]) + geno[,"Sex"])
  corrected <- rep(NA, length(geno[,tp]))
  notNA <- as.numeric(names(model$residuals))
  corrected[notNA] <- model$coefficients['(Intercept)'] + model$residuals
  return(corrected)
}

#############################################
#               Reading data                #
#############################################
setwd("D:/Edrive/Mouse/ClassicalPhenotypes/Complementation")

months <- rbind(c("Jan",31),c("Feb",28),c("Mar",31),c("Apr",30),c("May",31),c("Jun",30),
                c("Jul",31),c("Aug",31),c("Sep",30),c("Oct",31),c("Nov",30),c("Dec",31))

# AKR
AKR_MRI <- read.table("input/AKR_MRI.txt", sep='\t', header=TRUE, row.names=1)
AKR_MRI <- cbind(AKR_MRI, Age = NA)
AKR_Descr <- read.table("input/AKR_Descr.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", "x", "-"))

# NZO
NZO_Descr <- read.table("input/NZO_Descr.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", "x", "-"))
rownames(NZO_Descr) <-  gsub("777 ", "NZO-", rownames(NZO_Descr))

NZO_MRI <- read.table("input/NZO_MRI.txt", sep='\t', header=TRUE, row.names=1)
NZO_MRI <- cbind(NZO_MRI, Age = NA)

# TRPC
TRPC_Descr <- read.table("input/TRPC_Descr.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", "x", "-"))
rownames(TRPC_Descr) <-  gsub("TRPC", "TP", rownames(TRPC_Descr))

TRPC_MRI <- read.table("input/TRPC_MRI.txt", sep='\t', header=TRUE, row.names=1)
TRPC_MRI <- cbind(TRPC_MRI, Age = NA)

# BBS7
BBS7_Descr <- read.table("input/Bbs7_Descr.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", "x", "-"))
rownames(BBS7_Descr) <-  gsub("860BBS7-", "", rownames(BBS7_Descr))

BBS7_MRI <- read.csv("input/Bbs7_MRI.txt", sep='\t', header=TRUE, row.names=1)
BBS7_MRI <- cbind(BBS7_MRI, Age = NA)

#############################################
#             create MRI matrices           #
#############################################
fatlean_AKR <- createMRItable(AKR_MRI, AKR_Descr)
write.table(fatlean_AKR[[1]], "analysis/mri_fat_AKR.txt", sep = "\t")
write.table(fatlean_AKR[[2]], "analysis/mri_lean_AKR.txt", sep = "\t")

fatlean_NZO <- createMRItable(NZO_MRI, NZO_Descr)
write.table(fatlean_NZO[[1]], "analysis/mri_fat_NZO.txt", sep = "\t")
write.table(fatlean_NZO[[2]], "analysis/mri_lean_NZO.txt", sep = "\t")

fatlean_TRPC <- createMRItable(TRPC_MRI, TRPC_Descr)
write.table(fatlean_TRPC[[1]], "analysis/mri_fat_TRPC.txt", sep = "\t")
write.table(fatlean_TRPC[[2]], "analysis/mri_lean_TRPC.txt", sep = "\t")

fatlean_BBS7 <- createMRItable(BBS7_MRI, BBS7_Descr)
write.table(fatlean_BBS7[[1]], "analysis/mri_fat_BBS7.txt", sep = "\t")
write.table(fatlean_BBS7[[2]], "analysis/mri_lean_BBS7.txt", sep = "\t")

#############################################
#             filter individuals            #
#############################################
NZO_geno <- NZO_Descr[which(!is.na(NZO_Descr[,"C3GAB5"])),]
TRPC_geno <- TRPC_Descr[which(!is.na(TRPC_Descr[,"C3GAB5"])),]
AKR_geno <- AKR_Descr[which(!is.na(AKR_Descr[,"C3GAB5"])),]
BBS7_geno <- BBS7_Descr[which(!is.na(BBS7_Descr[,"C3GAB5"])),]

#############################################
#             add MRI day 70                #
#############################################
AKR_geno <- cbind(AKR_geno, fat70 = NA)
ii <- which(rownames(fatlean_AKR[[1]]) %in% rownames(AKR_geno))
AKR_geno[rownames(fatlean_AKR[[1]])[ii],"fat70"] <- fatlean_AKR[[1]][ii, "70"]

NZO_geno <- cbind(NZO_geno, fat70 = NA)
ii <- which(rownames(fatlean_NZO[[1]]) %in% rownames(NZO_geno))
NZO_geno[rownames(fatlean_NZO[[1]])[ii],"fat70"] <- fatlean_NZO[[1]][ii, "70"]

TRPC_geno <- cbind(TRPC_geno, fat70 = NA)
ii <- which(rownames(fatlean_TRPC[[1]]) %in% rownames(TRPC_geno))
TRPC_geno[rownames(fatlean_TRPC[[1]])[ii],"fat70"] <- fatlean_TRPC[[1]][ii, "70"]

BBS7_geno <- cbind(BBS7_geno, fat70 = NA)
ii <- which(rownames(fatlean_BBS7[[1]]) %in% rownames(BBS7_geno))
BBS7_geno[rownames(fatlean_BBS7[[1]])[ii],"fat70"] <- fatlean_BBS7[[1]][ii, "70"]

#############################################
#             recode genotypes              #
#############################################
gtypes <- c(193, 191, 189, 185, 177)
akrgname <- c("NZO", "B6N", "B6N", "BFMI", "AKR")
trpcgname <- c("NZO", "B6N", "B6N", "BFMI", "TRPC")

AKR.cg <- as.factor(unlist(lapply(lapply(strsplit(as.character(AKR_geno[,"C3GAB5"]), "/"),sort), function(x){
  paste0(akrgname[which(gtypes == x[1])],"|", akrgname[which(gtypes == x[2])])
})))

AKR_geno <- cbind(AKR_geno, Gtype = AKR.cg)

TRPC.cg <- as.factor(unlist(lapply(lapply(strsplit(as.character(TRPC_geno[,"C3GAB5"]), "/"),sort), function(x){
  paste0(trpcgname[which(gtypes == x[1])],"|", trpcgname[which(gtypes == x[2])])
})))

TRPC_geno <- cbind(TRPC_geno, Gtype = TRPC.cg)

NZO.cg <- as.factor(unlist(lapply(lapply(strsplit(as.character(NZO_geno[,"C3GAB5"]), "/"),sort), function(x){
  paste0(akrgname[which(gtypes == x[1])],"|", akrgname[which(gtypes == x[2])])
})))

NZO_geno <- cbind(NZO_geno, Gtype = NZO.cg)

#############################################
#        adjust phenotype and plot          #
#############################################

op <- par(mfrow=c(1,3))
boxplot(adjust.phenotype(AKR_geno, "63") ~ AKR_geno[,"Gtype"], main="F1 - AKR x (BFMI x B6N)", sub=paste0("Fat Day 63"))
boxplot(adjust.phenotype(NZO_geno, "63") ~ NZO_geno[,"Gtype"], main="F1 - NZO x (BFMI x B6N)", sub=paste0("Fat Day 63"))
boxplot(adjust.phenotype(TRPC_geno, "63") ~ TRPC_geno[,"Gtype"], main="F1 - TRPC x (BFMI x B6N)", sub=paste0("Fat Day 63"))

op <- par(mfrow=c(4,4))
for(x in colnames(NZO_geno)[13:16]) {
  boxplot(adjust.phenotype(NZO_geno, x) ~ NZO_geno[,"Gtype"], main="F1 - NZO x (BFMI x B6N)", xlab="", sub=paste0("Day ",x))
  cat("NZO:", anova(lm(adjust.phenotype(NZO_geno, x) ~ NZO_geno[,"Gtype"]))[[5]][1], "\n")
  boxplot(adjust.phenotype(TRPC_geno, x) ~ TRPC_geno[,"Gtype"], main="F1 - TRPC x (BFMI x B6N)", xlab="", sub=paste0("Day ",x))
  cat("TRPC:", anova(lm(adjust.phenotype(TRPC_geno, x) ~ TRPC_geno[,"Gtype"]))[[5]][1], "\n")
  boxplot(adjust.phenotype(AKR_geno, x) ~ AKR_geno[,"Gtype"], main="F1 - AKR x (BFMI x B6N)", xlab="", sub=paste0("Day ",x))
  cat("AKR:", anova(lm(adjust.phenotype(AKR_geno, x) ~ AKR_geno[,"Gtype"]))[[5]][1], "\n")
  boxplot(adjust.phenotype(BBS7_geno, x) ~ BBS7_geno[,"C3GAB5"], main="F1 - Bbs7 x (BFMI x B6N)", sub=paste0("Day ",x), xlab="", las=2)
  cat("Bbs7:", anova(lm(adjust.phenotype(BBS7_geno, x) ~ BBS7_geno[,"C3GAB5"]))[[5]][1], "\n")
  cat("Done",x,"\n")
  #line <- readline()
}

#############################################
#           use additional tests            #
#############################################

for(x in colnames(NZO_geno)[9:16]) {
  phe.adj <- adjust.phenotype(NZO_geno, x)
  iiB6N <- grepl("B6N", NZO_geno[,"Gtype"])
  iiBFMI <- grepl("BFMI", NZO_geno[,"Gtype"])
  cat("NZO t.test:", t.test(phe.adj[iiB6N], phe.adj[iiBFMI], alternative = "less")$p.value, "\n")  
}

for(x in colnames(AKR_geno)[9:16]) {
  phe.adj <- adjust.phenotype(AKR_geno, x)
  iiB6N <- grepl("B6N", AKR_geno[,"Gtype"])
  iiBFMI <- grepl("BFMI", AKR_geno[,"Gtype"])
  cat("AKR t.test:", t.test(phe.adj[iiB6N], phe.adj[iiBFMI], alternative = "less")$p.value, "\n")  
}

for(x in colnames(TRPC_geno)[9:16]) {
  phe.adj <- adjust.phenotype(TRPC_geno, x)
  iiB6N <- grepl("B6N", TRPC_geno[,"Gtype"])
  iiBFMI <- grepl("BFMI", TRPC_geno[,"Gtype"])
  cat("TRPC t.test:", t.test(phe.adj[iiB6N], phe.adj[iiBFMI], alternative = "less")$p.value, "\n")  
}

for(x in colnames(BBS7_geno)[9:16]) {
  phe.adj <- adjust.phenotype(BBS7_geno, x)
  iiB6N <- grepl("BBS7/+", BBS7_geno[,"C3GAB5"])
  iiBFMI <- grepl("BBS7/BFMI", BBS7_geno[,"C3GAB5"])
  cat("BBS7 t.test:", t.test(phe.adj[iiB6N], phe.adj[iiBFMI], alternative = "less")$p.value, "\n")  
}
