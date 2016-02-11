setwd("E:/Mouse/ClassicalPhenotypes/Complementation")

months <- rbind(c("Jan",31),c("Feb",28),c("Mar",31),c("Apr",30),c("May",31),c("Jun",30),
                c("Jul",31),c("Aug",31),c("Sep",30),c("Oct",31),c("Nov",30),c("Dec",31))


AKR_MRI <- read.table("input/AKR_MRI.txt", sep='\t', header=TRUE, row.names=1)
AKR_MRI <- rbind(AKR_MRI, read.table("input/AKR_MRI_Heise.txt", sep='\t', header=TRUE, row.names=1))
AKR_MRI <- cbind(AKR_MRI, Age = NA)

AKR_Descr <- read.table("input/AKR_Descr.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", "x", "-"))

source("D:\\Github\\HU-Berlin\\Mouse\\Muga\\Complementation\\createMRI.R")

fatlean_AKR <- createMRItable(AKR_MRI, AKR_Descr)
write.table(fatlean_AKR[[1]], "analysis/mri_fat_AKR.txt", sep = "\t")
write.table(fatlean_AKR[[2]], "analysis/mri_lean_AKR.txt", sep = "\t")

NZO_Descr <- read.table("input/NZO_Descr.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", "x", "-"))
NZO_MRI <- read.table("input/NZO_MRI.txt", sep='\t', header=TRUE, row.names=1)
NZO_MRI <- cbind(NZO_MRI, Age = NA)

fatlean_NZO <- createMRItable(NZO_MRI, NZO_Descr)
write.table(fatlean_NZO[[1]], "analysis/mri_fat_NZO.txt", sep = "\t")
write.table(fatlean_NZO[[2]], "analysis/mri_lean_NZO.txt", sep = "\t")

TRPC_Descr <- read.table("input/TRPC_Descr.txt", sep='\t', header=TRUE, row.names=1, check.names=FALSE, na.strings=c("NA", "", "x", "-"))

NZO_geno <- NZO_Descr[which(!is.na(NZO_Descr[,"C3GAB5"])),]
TRPC_geno <- TRPC_Descr[which(!is.na(TRPC_Descr[,"C3GAB5"])),]
AKR_geno <- AKR_Descr[which(!is.na(AKR_Descr[,"C3GAB5"])),]

AKR_Geno <- cbind(AKR_geno, fat70 = NA)
ii <- which(rownames(fatlean_AKR[[1]]) %in% rownames(AKR_Geno))
AKR_Geno[rownames(fatlean_AKR[[1]])[ii],"fat70"] <- fatlean_AKR[[1]][ii,"70"]

NZO_geno <- cbind(NZO_geno, fat70 = NA)
ii <- which(rownames(fatlean_NZO[[1]]) %in% rownames(NZO_geno))
NZO_geno[rownames(fatlean_NZO[[1]])[ii],"fat70"] <- fatlean_NZO[[1]][ii,"70"]


adjust.phenotype <- function(geno, tp = "21"){
  model <-  lm(geno[,tp] ~ geno[,"WG"] + geno[,"Sex"])
  corrected <- rep(NA, length(geno[,tp]))
  notNA <- as.numeric(names(model$residuals))
  corrected[notNA] <- model$coefficients['(Intercept)'] + model$residuals
  return(corrected)
}

op <- par(mfrow=c(1,2))
boxplot(adjust.phenotype(AKR_Geno, "fat70") ~ AKR_Geno[,"C3GAB5"], main="F1 - AKR x (BFMI x B6N)", sub=paste0("Fat Day 70"))
boxplot(adjust.phenotype(NZO_geno, "fat70") ~ NZO_geno[,"C3GAB5"], main="F1 - NZO x (BFMI x B6N)", sub=paste0("Fat Day 70"))

for(x in colnames(NZO_geno)[9:16]){
  op <- par(mfrow=c(1,3))
  boxplot(adjust.phenotype(NZO_geno, x) ~ NZO_geno[,"C3GAB5"], main="F1 - NZO x (BFMI x B6N)", sub=paste0("Day ",x))
  cat("NZO:", anova(lm(adjust.phenotype(NZO_geno, x) ~ NZO_geno[,"C3GAB5"]))[[5]][1], "\n")
  boxplot(adjust.phenotype(TRPC_geno, x) ~ TRPC_geno[,"C3GAB5"], main="F1 - TRPC x (BFMI x B6N)", sub=paste0("Day ",x))
  cat("TRPC:", anova(lm(adjust.phenotype(TRPC_geno, x) ~ TRPC_geno[,"C3GAB5"]))[[5]][1], "\n")
  boxplot(adjust.phenotype(AKR_geno, x) ~ AKR_geno[,"C3GAB5"], main="F1 - AKR x (BFMI x B6N)", sub=paste0("Day ",x))
  cat("AKR:", anova(lm(adjust.phenotype(AKR_geno, x) ~ AKR_geno[,"C3GAB5"]))[[5]][1], "\n")
  cat("Done",x,"\n")
  line <- readline()
}

