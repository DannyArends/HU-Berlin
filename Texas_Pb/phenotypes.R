setwd("D:/Ddrive/Github/HU-Berlin/Texas_Pb")
source("functions.R")

setwd("D:/Edrive/Mouse/Texas_Pb")
gts <- read.table("genotypes_F2_filtered_ordered.txt",sep="\t", header=TRUE, row.names=1)
phe <- read.table("input/F2_phenotypes.txt",sep="\t", header=TRUE, row.names=1, na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
rownames(phe) <- gsub("-", ".",rownames(phe))
phe <- phe[colnames(gts),]

# rename the column names, to be more logical and computer use-a-ble (no dots/spaces/etc)
colnames(phe) <- c("sex", "BLL", "DOB", "motherID1", "motherID2", "fatherID", "breedingCage", "weightBeforeDiet", "weightInitialMRI", 
                   "fatInitial", "leanInitial", "totalWaterInitial", "fatMassInitial", "leanMassInitial", "weightFinalMRI", "fatFinal", 
                   "leanFinal", "totalWaterFinal", "fatMassFinal", "leanMassFinal", "urineVolume", "waterConsumed", "waterConsumedDayMouse", 
                   "DSD", "ageAtDosing", "ageAtDosingW", "ageAtUC", "ageAtUCW", "fatGain")

phe <- cbind(phe, weightGain = phe[, "weightFinalMRI"] - phe[, "weightInitialMRI"])


phenotypes <- c("BLL", "weightBeforeDiet", "weightInitialMRI", "fatInitial", "leanInitial", "totalWaterInitial", "fatMassInitial", 
                "leanMassInitial", "weightFinalMRI", "fatFinal", "leanFinal", "totalWaterFinal", "fatMassFinal", "leanMassFinal", 
                "urineVolume", "waterConsumed", "waterConsumedDayMouse", "ageAtDosing", "ageAtDosingW", "ageAtUC", "ageAtUCW", "fatGain", "weightGain")

pheM <- phe[phe[, "sex"] == "M",phenotypes]
pheF <- phe[phe[, "sex"] == "F",phenotypes]

op <- par(mar=c(12,4,2,2))
plot(c(1,length(phenotypes)), c(-1,1), main = "Correlations",xaxt='n',xlab="", ylab="Pearson's R", las=2, t = 'n')
abline(h = seq(-1,1,0.25), lty=2, col="gray")
abline(h = 0, lty=2, col="black")
points(cor(pheM,use="pair")["BLL",], col="blue", pch=19,cex=2)
points(cor(pheF,use="pair")["BLL",], col="hotpink", pch=18,cex=2)
axis(1, at = 1:length(phenotypes), phenotypes, las=2)

# Since sex is a major driver we use the sex-corrected BLL to look for outliers:
BLLadj <- residuals(lm(phe[,"BLL"] ~ phe[, "sex"]))
hist(BLLadj, breaks = 30)

# We NA mask measurements > 3 SD from the mean
outlierBLL <- which(BLLadj > mean(BLLadj) + 3 * sd(BLLadj) | BLLadj < mean(BLLadj) - 3 * sd(BLLadj))
length(outlierBLL)
phe[outlierBLL,"BLL"] <- NA

# Since sex is a major driver we use the sex-corrected urineVolume to look for outliers:
urineVolumeadj <- residuals(lm(phe[,"urineVolume"] ~ phe[, "sex"]))
hist(urineVolumeadj, breaks = 30)

# We NA mask measurements > 3 SD from the mean
outlierUV <- which(urineVolumeadj > mean(urineVolumeadj) + 3 * sd(urineVolumeadj) | urineVolumeadj < mean(urineVolumeadj) - 3 * sd(urineVolumeadj))
length(outlierUV)
phe[outlierUV,"urineVolume"] <- NA

# Since sex is a major driver we use the sex-corrected waterConsumed to look for outliers:
waterConsumedadj <- residuals(lm(phe[,"waterConsumed"] ~ phe[, "sex"]))
hist(waterConsumedadj, breaks = 30)

# We NA mask measurements > 3 SD from the mean
outlierWC <- which(waterConsumedadj > mean(waterConsumedadj) + 3 * sd(waterConsumedadj) | waterConsumedadj < mean(waterConsumedadj) - 3 * sd(waterConsumedadj))
length(outlierWC)
phe[outlierWC,"waterConsumed"] <- NA

write.table(phe, "F2_phenotypes_cleaned_matched.txt", sep="\t", quote=FALSE)

f1 <- read.table("input/F1_phenotypes.txt",sep="\t", header=TRUE, row.names=1, na.strings=c("","-", "na", "NA", "NaN", "X", "x"))

# rename the column names, to be more logical and computer use-a-ble (no dots/spaces/etc)
colnames(f1) <- c("sex", "BLL", "DOB", "motherID1", "fatherID", "weightBeforeDiet", "weightInitialMRI", 
                   "fatInitial", "leanInitial", "totalWaterInitial", "fatMassInitial", "leanMassInitial", "weightFinalMRI", "fatFinal", 
                   "leanFinal", "totalWaterFinal", "fatMassFinal", "leanMassFinal", "urineVolume", "waterConsumed", "waterConsumedDayMouse", 
                   "DSD", "ageAtDosing", "ageAtDosingW", "ageAtUC", "ageAtUCW", "fatGain")

f1 <- cbind(f1, weightGain = f1[, "weightFinalMRI"] - f1[, "weightInitialMRI"])

computeH2 <- function(f1, f2){
  env <- mean((f1 - mean(f1, na.rm=TRUE))^2,na.rm=TRUE)
  total <- mean((f2 - mean(f2, na.rm=TRUE))^2,na.rm=TRUE)
  gen = round(100 * ((total - env) / total),1)
  return(gen)
}

pheM <- phe[phe[, "sex"] == "M",phenotypes]
pheF <- phe[phe[, "sex"] == "F",phenotypes]
f1M <- f1[f1[, "sex"] == "M",phenotypes]
f1F <- f1[f1[, "sex"] == "F",phenotypes]

gen.M <- c()
for(x in 1:1000){
  toR2 <- sample(which(!is.na(pheM[, "BLL"])),1)
  toR1 <- sample(which(!is.na(f1M[, "BLL"])),1)
  gen.M <- c(gen.M, computeH2(f1M[-toR1, "BLL"], pheM[-toR2, "BLL"]))
}
gen.F <- c()
for(x in 1:1000){
  toR2 <- sample(which(!is.na(pheF[, "BLL"])),1)
  toR1 <- sample(which(!is.na(f1F[, "BLL"])),1)
  gen.F <- c(gen.F, computeH2(f1F[-toR1, "BLL"], pheF[-toR2, "BLL"]))
}
mean(gen.M);sd(gen.M)
mean(gen.F);sd(gen.F)

gen.M <- c()
for(x in 1:1000){
  toR2 <- sample(which(!is.na(pheM[, "urineVolume"])),1)
  toR1 <- sample(which(!is.na(f1M[, "urineVolume"])),1)
  gen.M <- c(gen.M, computeH2(f1M[-toR1, "urineVolume"], pheM[-toR2, "urineVolume"]))
}
gen.F <- c()
for(x in 1:1000){
  toR2 <- sample(which(!is.na(pheF[, "urineVolume"])),1)
  toR1 <- sample(which(!is.na(f1F[, "urineVolume"])),1)
  gen.F <- c(gen.F, computeH2(f1F[-toR1, "urineVolume"], pheF[-toR2, "urineVolume"]))
}
mean(gen.M);sd(gen.M)
mean(gen.F);sd(gen.F)


gen.M <- c()
for(x in 1:1000){
  toR2 <- sample(which(!is.na(pheM[, "waterConsumed"])),1)
  toR1 <- sample(which(!is.na(f1M[, "waterConsumed"])),1)
  gen.M <- c(gen.M, computeH2(f1M[-toR1, "waterConsumed"], pheM[-toR2, "waterConsumed"]))
}
gen.F <- c()
for(x in 1:1000){
  toR2 <- sample(which(!is.na(pheF[, "waterConsumed"])),1)
  toR1 <- sample(which(!is.na(f1F[, "waterConsumed"])),1)
  gen.F <- c(gen.F, computeH2(f1F[-toR1, "waterConsumed"], pheF[-toR2, "waterConsumed"]))
}
mean(gen.M);sd(gen.M)
mean(gen.F);sd(gen.F)

