setwd("~/Github/HU-Berlin/Texas_Pb")
source("functions.R")

setwd("~/Data/Texas_Pb")
map <- read.table("reblasted_map.txt",sep="\t", header=TRUE, row.names=1)
gts <- read.table("genotypes_F2_filtered_ordered.txt",sep="\t", header=TRUE, row.names=1)
phe <- read.table("F2_phenotypes_cleaned_matched.txt",sep="\t", header=TRUE, row.names=1, na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
rownames(phe) <- gsub("-", ".",rownames(phe))
phe <- phe[colnames(gts),]
map <- map[rownames(gts),]

heatmap(cor(apply(phe[,c(2,8:23, 25:ncol(phe))],2,as.numeric), use="pair"))

# Some info
chroms <- c(1:19, "X")

dim(map)
dim(gts)
dim(phe)
phenotypes <- c("BLL", "weightBeforeDiet", "weightInitialMRI", "fatInitial", "leanInitial", "totalWaterInitial", "fatMassInitial", 
                "leanMassInitial", "weightFinalMRI", "fatFinal", "leanFinal", "totalWaterFinal", "fatMassFinal", "leanMassFinal", 
                "urineVolume", "waterConsumed", "waterConsumedDayMouse", "ageAtDosing", "ageAtDosingW", "ageAtUC", "ageAtUCW", "fatGain", "weightGain")

pheM <- phe[phe[, "sex"] == "M", phenotypes]
pheF <- phe[phe[, "sex"] == "F", phenotypes]

op <- par(mar=c(12,4,2,2))
plot(c(1,length(phenotypes)), c(-1,1), main = "Correlations",xaxt='n',xlab="", ylab="Pearson's R", las=2, t = 'n')
abline(h = seq(-1,1,0.25), lty=2, col="gray")
abline(h = 0, lty=2, col="black")
points(cor(pheM,use="pair")["BLL",], col="blue", pch=19,cex=2)
points(cor(pheF,use="pair")["BLL",], col="hotpink", pch=18,cex=2)
axis(1, at = 1:length(phenotypes), phenotypes, las=2)

# Minimal model (include sex as covariate)
anova(lm(BLL ~ sex, data = phe))
# Standard covars
anova(lm(BLL ~ sex + breedingCage, data = phe))
anova(lm(BLL ~ sex + ageAtDosing, data = phe))
anova(lm(BLL ~ sex + ageAtUC, data = phe))
# H2O related
anova(lm(BLL ~ sex + urineVolume, data = phe))
anova(lm(BLL ~ sex + waterConsumed, data = phe))
anova(lm(BLL ~ sex + totalWaterInitial, data = phe))
# Weights & MRI
anova(lm(BLL ~ sex + weightBeforeDiet, data = phe))
anova(lm(BLL ~ sex + weightInitialMRI, data = phe))
anova(lm(BLL ~ sex + weightFinalMRI, data = phe))
anova(lm(BLL ~ sex + fatInitial, data = phe))
anova(lm(BLL ~ sex + leanInitial, data = phe))
anova(lm(BLL ~ sex + fatGain, data = phe))
anova(lm(BLL ~ sex + weightGain, data = phe))

# Building a minimal model for QTL mapping of BLL
anova(lm(BLL ~ sex + urineVolume + ageAtDosing + ageAtUC + waterConsumed + weightBeforeDiet + weightInitialMRI+ weightFinalMRI+ breedingCage + weightInitialMRI + fatInitial + leanInitial + totalWaterInitial + fatGain, data = phe))
anova(lm(BLL ~ sex + urineVolume + waterConsumed + weightInitialMRI + leanInitial + totalWaterInitial, data = phe))
anova(lm(BLL ~ sex + urineVolume + waterConsumed + leanInitial + totalWaterInitial, data = phe))
anova(lm(BLL ~ sex + urineVolume + waterConsumed + totalWaterInitial, data = phe))
anova(lm(BLL ~ sex + urineVolume + waterConsumed, data = phe)) # All factors highly significant

mHas <- which(lapply(apply(gts[,rownames(pheM)],1,table), length) > 1) # Not all markers seggregate in males
fHas <- which(lapply(apply(gts[,rownames(pheF)],1,table), length) > 1) # Not all markers seggregate in females

# Models (Raw, Adjusted using the minimal model, GT interaction with Sex)
models.raw <- apply(gts,1,function(gt){ return(lm(BLL ~ sex + gt, data = phe)) })
models.adjusted <- apply(gts,1,function(gt){ return(lm(BLL ~ sex + urineVolume + waterConsumed + gt, data = phe)) })
models.interaction <- apply(gts, 1, function(gt){ return(lm(BLL ~ sex + urineVolume + waterConsumed + gt + sex:gt, data = phe)) })

# Male / Female only models
models.males <- apply(gts[mHas, rownames(pheM)], 1, function(gt){ return(lm(BLL ~ urineVolume + waterConsumed + gt, data = pheM)) })
models.females <- apply(gts[fHas, rownames(pheF)], 1, function(gt){ return(lm(BLL ~ urineVolume + waterConsumed + gt, data = pheF)) })

# Some variables for plotting
chrLengths <- chrLength(map)
map <- addCumLength(map,chrLengths)
chrC <- chrCenters(map,chrLengths)
colz <- rep(c("orange", "gray"), length(chroms))
colM <- rep(c("deepskyblue2", "darkturquoise"), length(chroms))
colF <- rep(c("deeppink", "deeppink3"), length(chroms))
colz <- colz[1:length(chroms)]
colM <- colM[1:length(chroms)]
colF <- colF[1:length(chroms)]
names(colz) <- chroms
names(colM) <- chroms
names(colF) <- chroms

# Anovas and resulting LOD scores
lod.raw <- toLOD(models.raw)
lod.adjusted <- toLOD(models.adjusted)

lod.males <- rep(NA, length(lod.adjusted))
names(lod.males) <- names(lod.adjusted)
lod.females <- rep(NA, length(lod.adjusted))
names(lod.females) <- names(lod.adjusted)

lod.males[names(toLOD(models.males))] <- toLOD(models.males)
lod.females[names(toLOD(models.females))] <- toLOD(models.females)

lod.interaction <- toLOD(models.interaction, "sex:gt")

## Thresholds from simpleM: blocksize 240, n Independant tests: 493
thresholds <- -log10(c(0.1,0.05,0.01)/493)

# Make the plot (Raw model, Combined Model, Sex:GT Interaction)
op <- par(mfrow=c(3,1))
op <- par(mar=c(5, 4, 4, 2) + 0.1)

plot(c(0, max(map[,"cumPos"],na.rm=TRUE)), c(0, 15), main = "Raw: BLL ~ sex + gt", t='n', xlab="Chromosome", ylab="-log10(P)", xaxt='n', las=2)
for(chr in chroms){ 
  idx <- which(map[,1] == chr)
  points(map[idx,"cumPos"], lod.raw[idx], col=colz[chr], t = 'l')
}
abline(h = thresholds, col=c("red", "gold", "green"), lty=2)
legend("topright", legend = paste0("p = ",c(0.1,0.05,0.01)), col=c("red", "gold", "green"), lty=2)
axis(1, at = chrC, chroms)

plot(c(0, max(map[,"cumPos"],na.rm=TRUE)), c(0, 15), main = "Combined: BLL ~ sex + urineVolume + waterConsumed + gt", t='n', xlab="Chromosome", ylab="-log10(P)", xaxt='n', las=2)
for(chr in chroms){ 
  idx <- which(map[,1] == chr)
  points(map[idx,"cumPos"], lod.adjusted[idx], col=colz[chr], t = 'l')
}
abline(h = thresholds, col=c("red", "gold", "green"), lty=2)
legend("topright", legend = paste0("p = ",c(0.1,0.05,0.01)), col=c("red", "gold", "green"), lty=2)
axis(1, at = chrC, chroms)

plot(c(0, max(map[,"cumPos"],na.rm=TRUE)), c(0, 15), main = "Interaction: BLL ~ sex + urineVolume + waterConsumed + gt + gt:sex", t='n', xlab="Chromosome", ylab="-log10(P)", xaxt='n', las=2)
for(chr in chroms){ 
  idx <- which(map[,1] == chr)
  points(map[idx,"cumPos"], lod.interaction[idx], col=colz[chr], t = 'l')
}
abline(h = thresholds, col=c("red", "gold", "green"), lty=2)
legend("topright", legend = paste0("p = ",c(0.1,0.05,0.01)), col=c("red", "gold", "green"), lty=2)
axis(1, at = chrC, chroms)

# Make the plot (Combined, Males, and Females)
op <- par(mfrow=c(3,1))
op <- par(mar=c(3, 4, 4, 2) + 0.1)
plot(c(0, max(map[,"cumPos"],na.rm=TRUE)), c(0, 12), main = "F2 population: BLL ~ sex + urineVolume + waterConsumed + gt", t='n', xlab="Chromosome", ylab="-log10(P)", xaxt='n', las=2)
for(chr in chroms){ 
  idx <- which(map[,1] == chr)
  points(map[idx,"cumPos"], lod.adjusted[idx], col=colz[chr], t = 'l')
}
abline(h = thresholds, col=c("red", "gold", "green"), lty=2)
legend("topright", legend = paste0("p = ",c(0.1,0.05,0.01)), col=c("red", "gold", "green"), lty=2)
axis(1, at = chrC, chroms)

plot(c(0, max(map[,"cumPos"],na.rm=TRUE)), c(0, 10), main = "Males: BLL ~ urineVolume + waterConsumed + gt", t='n', xlab="Chromosome", ylab="-log10(P)", xaxt='n', las=2)
for(chr in chroms){ 
  idx <- which(map[,1] == chr)
  points(map[idx,"cumPos"], lod.males[idx], col=colM[chr], t = 'l',lwd=2)
}
abline(h = thresholds, col=c("red", "gold", "green"), lty=2)
legend("topright", legend = paste0("p = ",c(0.1,0.05,0.01)), col=c("red", "gold", "green"), lty=2)
axis(1, at = chrC, chroms)

plot(c(0, max(map[,"cumPos"],na.rm=TRUE)), c(0, 10), main = "Females: BLL ~ urineVolume + waterConsumed + gt", t='n', xlab="Chromosome", ylab="-log10(P)", xaxt='n', las=2)
for(chr in chroms){ 
  idx <- which(map[,1] == chr)
  points(map[idx,"cumPos"], lod.females[idx], col=colF[chr], t = 'l',lwd=2)
}
abline(h = thresholds, col=c("red", "gold", "green"), lty=2)
legend("topright", legend = paste0("p = ",c(0.1,0.05,0.01)), col=c("red", "gold", "green"), lty=2)
axis(1, at = chrC, chroms)

summ.adj <- cbind(map[names(lod.adjusted),], lod.adjusted)
summ.adj[which(summ.adj[,1] == 1),]
summ.adj[which(summ.adj[,1] == 7),]

chr1Top <- "gUNC2018619"
chr7Top <- "gUNC13158239"

gts.num <- read.table("genotypes_F2_filtered_ordered_numeric.txt", sep = "\t")

BLLadj <- residuals(lm(BLL ~ sex + urineVolume + waterConsumed, data = phe)) + mean(phe[,"BLL"], na.rm=TRUE)

phe <- cbind(phe[,1:2], BLLadj = round(BLLadj[rownames(phe)],1), phe[,3:ncol(phe)])
write.table(phe, "phe_BLLadjusted.txt",sep="\t", quote=FALSE)

op <- par(mfrow=c(1,2))
plot(c(0.4, 3.6), c(0,100), t = 'n', xaxt='n', las=2, ylab = "Adjusted blood lead level (BLL)", xlab="Genotype", main=paste0("BLL - Combined - ",chr1Top))
boxplot(BLLadj ~ as.character(gts[chr1Top,names(BLLadj)]), notch=TRUE, yaxt='n', add=TRUE, col=c("gray25", "gray50", "gray75"))
plot(c(0.4, 3.6), c(0,100), t = 'n', xaxt='n', las=2, ylab = "Adjusted blood lead level (BLL)", xlab="Genotype", main=paste0("BLL - Combined - ",chr7Top))
boxplot(BLLadj ~ as.character(gts[chr7Top,names(BLLadj)]), notch=TRUE, yaxt='n', add=TRUE, col=c("gray25", "gray50", "gray75"))

summ.males <- cbind(map[names(lod.males),], lod.males)
summ.males[which(summ.males[,1] == 6),]
summ.males[which(summ.males[,1] == 7),]

summ.females <- cbind(map[names(lod.females),], lod.females)
summ.females[which(summ.females[,1] == 1),]
summ.females[which(summ.females[,1] == 7),]

### Pairscan
pairscan <- matrix(NA, nrow(gts), nrow(gts), dimnames = list(rownames(gts), rownames(gts)))

#i <- 2
#for(m1 in rownames(gts)){
#  for(m2 in rownames(gts)[i:nrow(gts)]){
#    gt1 <- as.character(gts[m1,])
#    gt2 <- as.character(gts[m2,])
#    pairscan[m1,m2] <- anova(lm(BLL ~ sex + urineVolume + waterConsumed + gt1 + gt2 + gt1:gt2, data = phe))["gt1:gt2", "Pr(>F)"]
#  }
#  cat("Finished, with marker",i,"\n")
#  i <- i + 1
#}
#write.table(pairscan, "results.pairscan.txt", sep = "\t", quote=FALSE)


setwd("D:/Edrive/Mouse/Texas_Pb")
map <- read.table("reblasted_map.txt",sep="\t", header=TRUE, row.names=1)
gts <- read.table("genotypes_F2_filtered_ordered.txt",sep="\t", header=TRUE, row.names=1)
phe <- read.table("F2_phenotypes_cleaned_matched.txt",sep="\t", header=TRUE, row.names=1, na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
rownames(phe) <- gsub("-", ".",rownames(phe))
phe <- phe[colnames(gts),]
map <- map[rownames(gts),]
pairscan <- read.table("results.pairscan.txt", sep = "\t")
pairscan <- pairscan[,-1]


threshold <- -log10(c(0.1, 0.05, 0.01)/(493 * 493))

pairscanLODs <- as.matrix(-log10(pairscan))
image(1:ncol(pairscanLODs), 1:nrow(pairscanLODs), t(pairscanLODs),xlab="Chromosome", ylab="Chromosome", xaxt='n', yaxt='n', breaks = c(0, 3, 6, 7,10), col=c("gray", "yellow", "orange", "red"))
tbls <- table(map[rownames(pairscanLODs),1])
chrs <- c()
for(chr in chroms){
  chrs <- c(chrs, chr)
  abline(h = sum(tbls[chrs]), lty=2)
  abline(v = sum(tbls[chrs]), lty=2)
}

m1 <- rownames(map[which(apply(pairscanLODs,1, function(x){any(x > 5,na.rm=TRUE)})),])
m18 <- rownames(map[which(apply(pairscanLODs,2, function(x){any(x > 5,na.rm=TRUE)})),])

pairscanLODs[m1,m18]

BLLadj <- mean(phe[, "BLL"],na.rm=TRUE) + residuals(lm(BLL ~ sex + urineVolume + waterConsumed, data = phe))
top1 <- "gJAX00272462"
top18 <- "gUNC29787730"
topX <- "SXX202459414"
chr1Top <- "backupUNC010628673"
chr7Top <- "gUNC13158239"

op <- par(mfrow=c(1,2))

plot(c(0.4, 3.6), c(0,100), t = 'n', xaxt='n', las=2, ylab = "Adjusted blood lead level (BLL)", xlab="Genotype backupUNC010628673", main=paste0("BLL - Combined - Chromosome 1"),cex.axis=1.5,cex.lab=1.5)
boxplot(BLLadj ~ as.character(gts[chr1Top,names(BLLadj)]), notch=TRUE, yaxt='n', add=TRUE, col=c("gray25", "gray50", "gray75"),cex.axis=1.5,cex.lab=1.5)
boxplot(BLLadj ~ as.character(gts[top18,names(BLLadj)]), notch=TRUE, yaxt='n', add=TRUE, col=c("gray25", "gray50", "gray75"),cex.axis=1.5,cex.lab=1.5)

plot(c(0.4, 3.6), c(0,100), t = 'n', xaxt='n', las=2, ylab = "Adjusted blood lead level (BLL)", xlab="Genotype gUNC13158239", main=paste0("BLL - Combined - Chromosome 7"),cex.axis=1.5,cex.lab=1.5)
boxplot(BLLadj ~ as.character(gts[chr7Top,names(BLLadj)]), notch=TRUE, yaxt='n', add=TRUE, col=c("gray25", "gray50", "gray75"),cex.axis=1.5,cex.lab=1.5)

op <- par(mfrow=c(1,1))

fem <- rownames(phe[which(phe[colnames(gts), "sex"] == "F"),])
males <- rownames(phe[which(phe[colnames(gts), "sex"] == "M"),])

layout(matrix(c(1, 1, 2, 3), 2, 2, byrow = FALSE))

plot(c(0.75,3.25), c(20, 70), t = 'n',xaxt='n', xlab="Genotype gJAX00272462 (chr1)", ylab = "Adjusted blood lead level (BLL)", xaxs = 'i', yaxs = 'i', main = "Interaction Chr1 and Chr18 (Combined)",cex.axis=1.5,cex.lab=1.2, las=2)
abline(h = seq(0,100,5), lty=2, col = "gray")
#abline(v = 1:3, lty=2, col = "gray")
cnt <- 1
for(gt18 in c("AA", "AG", "GG")){
  means <- c()
  at <- 1 + (cnt-2) / 20
  for(gt1 in c("AA", "AG", "GG")){
    ind <- colnames(gts)[which(gts[top1,] == gt1 & gts[top18,] == gt18)]
    mymean <- median(BLLadj[ind],na.rm=TRUE)
    myse <- sd(BLLadj[ind],na.rm=TRUE)/sqrt(length(ind))
    means <- c(means, mymean)
    points(c(at,at), c(mymean-myse, mymean+myse), col=cnt, t = 'l', lwd=2)
    points(c(at-0.02,at+0.015), c(mymean+myse, mymean+myse), col=cnt, t = 'l', lwd=1)
    points(c(at-0.02,at+0.015), c(mymean-myse, mymean-myse), col=cnt, t = 'l', lwd=1)
    at <- at + 1
    cat(gt1, gt18, "n = ",length(ind), mymean, "\n")
  }
  points(1:3 + (cnt-2) / 20, means, t = 'l', col = cnt, lwd=2)
  cnt <- cnt + 1
}
legend("topright", c("AA", "AG", "GG"), lty=1, title  = "Genotype gUNC29787730 (chr18)", col = 1:3, bg="white", ncol = 3)
axis(1, at = c(1,2,3), c("AA", "AG", "GG"),cex.axis=1.5,cex.lab=1.5)
box()

#op <- par(mfrow = c(1,2))

plot(c(0.75,3.25), c(20, 70), t = 'n',xaxt='n', xlab="Genotype gJAX00272462 (chr1)", ylab = "Adjusted blood lead level (BLL)", xaxs = 'i', yaxs = 'i', main = "Interaction Chr1 and Chr18 (Females)",cex.axis=1.2,cex.lab=1, las=2)
abline(h = seq(0,100,5), lty=2, col = "gray")
#abline(v = 1:3, lty=2, col = "gray")
cnt <- 1
for(gt18 in c("AA", "AG", "GG")){
  means <- c()
  at <- 1 + (cnt-2) / 20
  for(gt1 in c("AA", "AG", "GG")){
    ind <- colnames(gts)[which(gts[top1,] == gt1 & gts[top18,] == gt18)]
    ind <- ind[ind %in% fem]
    mymean <- median(BLLadj[ind],na.rm=TRUE)
    myse <- sd(BLLadj[ind],na.rm=TRUE)/sqrt(length(ind))
    means <- c(means, mymean)
    points(c(at,at), c(mymean-myse, mymean+myse), col=cnt, t = 'l', lwd=2)
    points(c(at-0.02,at+0.015), c(mymean+myse, mymean+myse), col=cnt, t = 'l', lwd=1)
    points(c(at-0.02,at+0.015), c(mymean-myse, mymean-myse), col=cnt, t = 'l', lwd=1)
    at <- at + 1
    cat(gt1, gt18, "n = ",length(ind), mymean, "\n")
  }
  points(1:3 + (cnt-2) / 20, means, t = 'l', col = cnt, lwd=2)
  cnt <- cnt + 1
}
legend("topright", c("AA", "AG", "GG"), lty=1, title  = "Genotype gUNC29787730 (chr18)", col = 1:3, bg="white", ncol = 3)
axis(1, at = c(1,2,3), c("AA", "AG", "GG"),cex.axis=1.5,cex.lab=1.5)
box()


plot(c(0.75,3.25), c(20, 70), t = 'n',xaxt='n', xlab="Genotype gJAX00272462 (chr1)", ylab = "Adjusted blood lead level (BLL)", xaxs = 'i', yaxs = 'i', main = "Interaction Chr1 and Chr18 (Males)",cex.axis=1.2,cex.lab=1, las=2)
abline(h = seq(0,100,5), lty=2, col = "gray")
#abline(v = 1:3, lty=2, col = "gray")
cnt <- 1
for(gt18 in c("AA", "AG", "GG")){
  means <- c()
  at <- 1 + (cnt-2) / 20
  for(gt1 in c("AA", "AG", "GG")){
    ind <- colnames(gts)[which(gts[top1,] == gt1 & gts[top18,] == gt18)]
    ind <- ind[ind %in% males]
    mymean <- median(BLLadj[ind],na.rm=TRUE)
    myse <- sd(BLLadj[ind],na.rm=TRUE)/sqrt(length(ind))
    means <- c(means, mymean)
    points(c(at,at), c(mymean-myse, mymean+myse), col=cnt, t = 'l', lwd=2)
    points(c(at-0.02,at+0.015), c(mymean+myse, mymean+myse), col=cnt, t = 'l', lwd=1)
    points(c(at-0.02,at+0.015), c(mymean-myse, mymean-myse), col=cnt, t = 'l', lwd=1)
    at <- at + 1
    cat(gt1, gt18, "n = ",length(ind), mymean, "\n")
  }
  points(1:3 + (cnt-2) / 20, means, t = 'l', col = cnt, lwd=2)
  cnt <- cnt + 1
}
legend("topright", c("AA", "AG", "GG"), lty=1, title  = "Genotype gUNC29787730 (chr18)", col = 1:3, bg="white", ncol = 3)
axis(1, at = c(1,2,3), c("AA", "AG", "GG"),cex.axis=1.5,cex.lab=1.5)
box()


xAA <- colnames(gts)[which(gts[topX,] == "AA")]
xAT <- colnames(gts)[which(gts[topX,] == "AT")]
xTT <- colnames(gts)[which(gts[topX,] == "TT")]

op <- par(mfrow = c(1,3))

plot(c(0.75,3.25), c(20, 70), t = 'n',xaxt='n', xlab="Genotype chromosome 1", ylab = "Adjusted blood lead level (BLL)", xaxs = 'i', yaxs = 'i', main = "Interaction Chr1 and Chr18 (Chr X=AA)",cex.axis=1.5,cex.lab=1, las=2)
abline(h = seq(0,100,5), lty=2, col = "gray")
#abline(v = 1:3, lty=2, col = "gray")
cnt <- 1
for(gt18 in c("AA", "AG", "GG")){
  means <- c()
  at <- 1 + (cnt-2) / 20
  for(gt1 in c("AA", "AG", "GG")){
    ind <- colnames(gts)[which(gts[top1,] == gt1 & gts[top18,] == gt18)]
    ind <- ind[ind %in% xAA]
    mymean <- median(BLLadj[ind],na.rm=TRUE)
    myse <- sd(BLLadj[ind],na.rm=TRUE)/sqrt(length(ind))
    means <- c(means, mymean)
    points(c(at,at), c(mymean-myse, mymean+myse), col=cnt, t = 'l', lwd=2)
    points(c(at-0.02,at+0.015), c(mymean+myse, mymean+myse), col=cnt, t = 'l', lwd=1)
    points(c(at-0.02,at+0.015), c(mymean-myse, mymean-myse), col=cnt, t = 'l', lwd=1)
    at <- at + 1
    cat(gt1," ", gt18,"n = ",length(ind), "\n")
  }
  points(1:3 + (cnt-2) / 20, means, t = 'l', col = cnt, lwd=2)
  cnt <- cnt + 1
}
legend("topleft", c("AA", "AG", "GG"), lty=1, title  = "Chromosome 18", col = 1:3, bg="white")
axis(1, at = c(1,2,3), c("AA", "AG", "GG"),cex.axis=1.5,cex.lab=1.5)
box()


plot(c(0.75,3.25), c(20, 70), t = 'n',xaxt='n', xlab="Genotype chromosome 1", ylab = "Adjusted blood lead level (BLL)", xaxs = 'i', yaxs = 'i', main = "Interaction Chr1 and Chr18 (Chr X = AT)",cex.axis=1.5,cex.lab=1.5, las=2)
abline(h = seq(0,100,5), lty=2, col = "gray")
#abline(v = 1:3, lty=2, col = "gray")
cnt <- 1
for(gt18 in c("AA", "AG", "GG")){
  means <- c()
  at <- 1 + (cnt-2) / 20
  for(gt1 in c("AA", "AG", "GG")){
    ind <- colnames(gts)[which(gts[top1,] == gt1 & gts[top18,] == gt18)]
    ind <- ind[ind %in% xAT]
    mymean <- median(BLLadj[ind],na.rm=TRUE)
    myse <- sd(BLLadj[ind],na.rm=TRUE)/sqrt(length(ind))
    means <- c(means, mymean)
    points(c(at,at), c(mymean-myse, mymean+myse), col=cnt, t = 'l', lwd=2)
    points(c(at-0.02,at+0.015), c(mymean+myse, mymean+myse), col=cnt, t = 'l', lwd=1)
    points(c(at-0.02,at+0.015), c(mymean-myse, mymean-myse), col=cnt, t = 'l', lwd=1)
    at <- at + 1
    cat("n = ",length(ind), "\n")
  }
  points(1:3 + (cnt-2) / 20, means, t = 'l', col = cnt, lwd=2)
  cnt <- cnt + 1
}
legend("topleft", c("AA", "AG", "GG"), lty=1, title  = "Chromosome 18", col = 1:3, bg="white")
axis(1, at = c(1,2,3), c("AA", "AG", "GG"),cex.axis=1.5,cex.lab=1.5)
box()

plot(c(0.75,3.25), c(20, 70), t = 'n',xaxt='n', xlab="Genotype chromosome 1", ylab = "Adjusted blood lead level (BLL)", xaxs = 'i', yaxs = 'i', main = "Interaction Chr1 and Chr18 (Chr X = TT)",cex.axis=1.5,cex.lab=1.5, las=2)
abline(h = seq(0,100,5), lty=2, col = "gray")
#abline(v = 1:3, lty=2, col = "gray")
cnt <- 1
for(gt18 in c("AA", "AG", "GG")){
  means <- c()
  at <- 1 + (cnt-2) / 20
  for(gt1 in c("AA", "AG", "GG")){
    ind <- colnames(gts)[which(gts[top1,] == gt1 & gts[top18,] == gt18)]
    ind <- ind[ind %in% xTT]
    mymean <- median(BLLadj[ind],na.rm=TRUE)
    myse <- sd(BLLadj[ind],na.rm=TRUE)/sqrt(length(ind))
    means <- c(means, mymean)
    points(c(at,at), c(mymean-myse, mymean+myse), col=cnt, t = 'l', lwd=2)
    points(c(at-0.02,at+0.015), c(mymean+myse, mymean+myse), col=cnt, t = 'l', lwd=1)
    points(c(at-0.02,at+0.015), c(mymean-myse, mymean-myse), col=cnt, t = 'l', lwd=1)
    at <- at + 1
    cat("n = ",length(ind), "\n")
  }
  points(1:3 + (cnt-2) / 20, means, t = 'l', col = cnt, lwd=2)
  cnt <- cnt + 1
}
legend("topleft", c("AA", "AG", "GG"), lty=1, title  = "Chromosome 18", col = 1:3, bg="white")
axis(1, at = c(1,2,3), c("AA", "AG", "GG"),cex.axis=1.5,cex.lab=1.5)
box()
