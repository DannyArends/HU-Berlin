# Data Quality control
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014
library(lme4)

source("D:/Github/HU-Berlin/Mouse/Muga/ATB_paper/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes
 
missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("33310233", "6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes  <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                         # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes

GinP <- which(colnames(genotypes) %in% rownames(phenotypes))
selected <- colnames(genotypes)[GinP]

write.table(cbind(rownames(genotypes), genotypes[,selected]),"D:/Papers/megaMuga/Additional files/Supplement2 - Genotypes.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(cbind(rownames(phenotypes), phenotypes)[selected,],"D:/Papers/megaMuga/Additional files/Supplement3 - Phenotypes.txt",sep="\t",quote=FALSE,row.names=FALSE)

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]                                                 # The F2 individuals

onegenotype <- which(lapply(apply(genotypes[,F2], 1, table), length) < 2)                                       # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- (genotypes[-onegenotype, F2])                                                                    # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypes), "markers\n")                                                                  # == Left with 21260 markers    // Left with 11581 markers

colorz <- as.numeric(as.factor(unlist(genotypes["UNC5048297", F2])))
growth <- phenotypes[F2, c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")]

fat <- phenotypes[F2, c("mri42d_fat", "mri56d_fat", "mri70d_fat")]
lean <- phenotypes[F2, c("mri42d_lean", "mri56d_lean", "mri70d_lean")]
fatlean <- (fat / (fat+lean))

vater <- as.factor(phenotypes[F2, c("Vater")]) ; mutter <- as.factor(phenotypes[F2, c("Mutter")]) ; wsize <- as.numeric(phenotypes[F2, c("WG2")]) ; wlabel <- as.factor(phenotypes[F2, c("W.Label")]) ; season <- as.factor(phenotypes[F2, c("Season")])
fam <- as.factor(paste0(mutter, vater))
#op <- par(mfrow=c(3,2))
#boxplot(growth[,"d70"] ~ colorz, main="without")

#nVgrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ vater ); return(aa$coefficients["(Intercept)"] + aa$residuals) })
#boxplot(nVgrowth[,"d70"] ~ colorz, main="Vater")
#nWLgrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wlabel); return(aa$coefficients["(Intercept)"] + aa$residuals) })
#boxplot(nWLgrowth[,"d70"] ~ colorz, main="wLabel")
#nSgrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ season); return(aa$coefficients["(Intercept)"] + aa$residuals) })
#boxplot(nSgrowth[,"d70"] ~ colorz, main="Season")
#nWSgrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wsize); return(aa$coefficients["(Intercept)"] + aa$residuals) })
#boxplot(nWSgrowth[,"d70"] ~ colorz, main="wSize")

#nWSFgrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ vater + wsize); return(aa$coefficients["(Intercept)"] + aa$residuals) })
#boxplot(nWSFgrowth[,"d70"] ~ colorz, main="father + wSize")


boxplot(growth)

boxplot(growth[colorz == 1,], col=rgb(1,0,0,0.5), main = "Raw data")
boxplot(growth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(growth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(growth[colorz == 1,]) - colMeans(growth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(growth[colorz == 1,]) - colMeans(growth[colorz == 3,]),d=2), "\n")

## Family (Fixed)

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ fam ); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Family")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Family (Random)

ngrowth <- apply(growth, 2, function(x){ aa <- lmer(as.numeric(x) ~ (1|fam) ); return(residuals(aa)) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Family")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Litter size

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wsize ); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Litter size")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Family + Litter size

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ fam + wsize ); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Family + Litter size")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Family per Geno

boxplot(cBfmi <- apply(growth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(1,0,0,0.5), main = "~ Family(Geno)")
boxplot(cHetro <- apply(growth[colorz == 2,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 2]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(cB6n <- apply(growth[colorz == 3,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 3]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(cHetro),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(cB6n),d=2), "\n")

## Family per Geno (merged B6n and Hetro)
cBfmi <- apply(growth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1]); return(aa$coefficients["(Intercept)"] + aa$residuals) })
merged <- apply(growth[colorz != 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz != 1]); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(cBfmi, col=rgb(1,0,0,0.5), main = "~ Family(Geno)")
boxplot(merged[colorz[colorz != 1] == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(merged[colorz[colorz != 1] == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(merged[colorz[colorz != 1] == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(merged[colorz[colorz != 1] == 3,]),d=2), "\n")

## Family and Litter size per Geno

boxplot(cBfmi <- apply(growth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1] + wsize[colorz==1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(1,0,0,0.5), main = "~ Family(Geno) + Litter size(Geno)")
boxplot(cHetro <- apply(growth[colorz == 2,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 2]+ wsize[colorz==2]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(cB6n <- apply(growth[colorz == 3,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 3]+ wsize[colorz==3]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(cHetro),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(cB6n),d=2), "\n")

## Family and Litter size per Geno (merged B6n and Hetro)
cBfmi <- apply(growth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1]+ wsize[colorz==1]); return(aa$coefficients["(Intercept)"] + aa$residuals) })
merged <- apply(growth[colorz != 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz != 1]+ wsize[colorz!=1]); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(cBfmi, col=rgb(1,0,0,0.5), main = "~ Family(Geno) + Litter size(Geno)")
boxplot(merged[colorz[colorz != 1] == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(merged[colorz[colorz != 1] == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(merged[colorz[colorz != 1] == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(merged[colorz[colorz != 1] == 3,]),d=2), "\n")

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
qtl42FatO  <- read.table("Analysis/qtls_mri42d_fat.txt",   sep="\t")
qtl42Fat   <- read.table("Analysis/kasp/QTL_mri42d_fat_kasp.txt",   sep="\t")
qtl42Fat   <- rbind(qtl42Fat, qtl42FatO)

setwd("E:/Mouse/DNA/MegaMuga/")

kasp   <- read.table("Analysis/kasp+licor.txt", sep="\t", check.names=FALSE, colClasses="character", skip=2, header=TRUE, na.strings=c("","-", "?","NA", "CG?", "CT?", "CT?", "AG?", "KL?"))

map <- t(kasp[c(1,3), -c(1:6)])
colnames(map) <- c("Chr", "Mb_NCBI38")

kasp <- kasp[-c(1:3), -c(2:6)]
rownames(kasp) <- kasp[,1]
kasp <- kasp[,-1]

goodM <- colnames(kasp)[ which(apply(kasp,2,function(x){return(length(which(is.na(x)))) }) < 400) ] # colnames of good markers

genotypes <- t(kasp[,goodM])   # Individuals in the columns, SNPs in the rows
map <- map[goodM,]        # SNPs in the rows, Chromosome and Position (NCBI38)

mapO <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
map <- rbind(map,mapO[,1:2])

ii <- rownames(qtl42Fat[qtl42Fat[,"marker"] > 1,])
submap <- map[ii, ]
submap <- submap[submap[,"Chr"] == 3,]                                    #  & 
submap <- submap[rownames(submap) %in% rownames(qtl42Fat),]


submap <- submap[as.numeric(as.character(submap[,"Mb_NCBI38"])) > 32500000 &  as.numeric(as.character(submap[,"Mb_NCBI38"])) < 40000000,]
submap <- submap[order(as.numeric(as.character(submap[,2]))),]

colz <- rep("black", nrow(submap))
colz[grep("KM", rownames(submap))] <- "red"

## Current model *No family effect*
op <- par(mfrow=c(2, 2), cex.main=1.3, cex.axis=1.2,cex.lab=1.3)

plot(c(32500000, 40000000), y = c(0, 65), t = 'n', ylab = "-log10(pvalue)", xlab = "Chromosome 3: 32.5 - 40 Mb", xaxt='n', las = 2, main= "a) QTL profile Fat mass (week 6)")
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = qtl42Fat[rownames(submap), "marker"], t = 'l', col = "red",lwd=3)
points(x = as.numeric(as.character(submap[,"Mb_NCBI38"])), y = rep(-1.3, nrow(submap)), pch="|", col = colz, cex=0.5)
axis(1, at=seq(32500000, 40000000, 2500000), seq(32500000, 40000000, 2500000) / 1000000)

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wsize + wlabel + season); return(aa$coefficients["(Intercept)"] + aa$residuals) })
colnames(ngrowth) <- seq(21,70,7)
plot(c(20,72),c(10,65), t = 'n', xaxt='n', main = "b) Body mass", ylab = "Body mass (g)", xlab="Age (weeks)",las=2, xaxt='n')
rect(0,0,((28-21)/2) + 21,100, col=rgb(0,0,1,0.3), border=FALSE)
rect(((42-35)/2) + 35,0,100,100, col=rgb(0,1,0,0.3), border=FALSE)
boxplot(ngrowth[colorz == 1,], at = seq(21,70,7)+1.1, col=rgb(1,0.4,0,0.5), pars=list(boxwex=.8), xaxt='n', notch=TRUE,add = TRUE, yaxt='n', xaxt='n')
boxplot(ngrowth[colorz == 2,], at = seq(21,70,7), col=rgb(0.5,0.5,0.5,1), add=TRUE, pars=list(boxwex=.8), notch=TRUE, yaxt='n', xaxt='n')
boxplot(ngrowth[colorz == 3,], at = seq(21,70,7)-1.1, col=rgb(0,0,1,0.5), add=TRUE, pars=list(boxwex=.8),xaxt='n', notch=TRUE, yaxt='n', xaxt='n')
points(x = seq(21,70,7), y = apply(ngrowth[colorz == 2,],2,median),t='l', col=rgb(0.5,0.5,0.5,1),lwd=2)
points(x = seq(21,70,7), y = apply(ngrowth[colorz == 3,],2,median),t='l', col=rgb(0,0,1,0.5),lwd=2)
points(x = seq(21,70,7), y = apply(ngrowth[colorz == 1,],2,median),t='l', col=rgb(1,0.4,0,0.5),lwd=2)
legend("topleft", c("BFMI","Heterozygous", "B6N"), fill = c(rgb(1,0.4,0,0.5),rgb(0.5,0.5,0.5,1),rgb(0,0,1,0.5)), title="UNC5048297 genotype",bg="white",cex=1.2)
axis(1, at= seq(21, 72, 7), seq(21, 72, 7) / 7)

nfat <- apply(fat, 2, function(x){ aa <- lm(as.numeric(x) ~ wsize + wlabel + season); bb <- rep(NA,length(x));bb[as.numeric(names(aa$residuals))] <- aa$residuals; bb <- bb + aa$coefficients["(Intercept)"]; return(bb) })
colnames(nfat) <- c(42, 56, 70)
plot(c(40,72),c(0,25), t = 'n', xaxt='n', main = "c) Fat mass", ylab = "Fat mass (g)", xlab="Age (weeks)",las=2, xaxt='n')
boxplot(nfat[colorz == 1,], at = c(42, 56, 70)+1.1, col=rgb(1,0.4,0,0.5), pars=list(boxwex=.8), xaxt='n', notch=TRUE,add = TRUE, yaxt='n', xaxt='n')
boxplot(nfat[colorz == 2,], at = c(42, 56, 70), col=rgb(0.5,0.5,0.5,1), add=TRUE, pars=list(boxwex=.8), notch=TRUE, yaxt='n', xaxt='n')
boxplot(nfat[colorz == 3,], at = c(42, 56, 70)-1.1, col=rgb(0,0,1,0.5), add=TRUE, pars=list(boxwex=.8),xaxt='n', notch=TRUE, yaxt='n', xaxt='n')
points(x = c(42, 56, 70), y = apply(nfat[colorz == 1,],2, median,na.rm=TRUE),t='l', col=rgb(1,0.4,0,0.5),lwd=2)
points(x = c(42, 56, 70), y = apply(nfat[colorz == 2,],2, median,na.rm=TRUE),t='l', col=rgb(0.5,0.5,0.5,1),lwd=2)
points(x = c(42, 56, 70), y = apply(nfat[colorz == 3,],2, median,na.rm=TRUE),t='l', col=rgb(0,0,1,0.5),lwd=2)
legend("topleft", c("BFMI","Heterozygous", "B6N"), fill = c(rgb(1,0.4,0,0.5),rgb(0.5,0.5,0.5,1),rgb(0,0,1,0.5)), title="UNC5048297 genotype",bg="white",cex=1.2)
axis(1, at= seq(14, 72, 14), seq(14, 72,14) / 7)

nfatlean <- apply(fatlean, 2, function(x){ aa <- lm(as.numeric(x) ~ wsize + wlabel + season); bb <- rep(NA,length(x));bb[as.numeric(names(aa$residuals))] <- aa$residuals; bb <- bb + aa$coefficients["(Intercept)"]; return(bb) })
colnames(nfatlean) <- c(42, 56, 70)
plot(c(40,72),c(0,0.75), t = 'n', xaxt='n', yaxt='n', main = "d) Fat percentage", ylab = "Fat / (Fat + Lean) (%)", xlab="Age (weeks)",las=2, xaxt='n')
boxplot(nfatlean[colorz == 1,], at = c(42, 56, 70)+1.1, col=rgb(1,0.4,0,0.5), pars=list(boxwex=.8), xaxt='n', notch=TRUE,add = TRUE, yaxt='n', xaxt='n')
boxplot(nfatlean[colorz == 2,], at = c(42, 56, 70), col=rgb(0.5,0.5,0.5,1), add=TRUE, pars=list(boxwex=.8), notch=TRUE, yaxt='n', xaxt='n')
boxplot(nfatlean[colorz == 3,], at = c(42, 56, 70)-1.1, col=rgb(0,0,1,0.5), add=TRUE, pars=list(boxwex=.8),xaxt='n', notch=TRUE, yaxt='n', xaxt='n')
points(x = c(42, 56, 70), y = apply(nfatlean[colorz == 1,],2, median,na.rm=TRUE),t='l', col=rgb(1,0.4,0,0.5),lwd=2)
points(x = c(42, 56, 70), y = apply(nfatlean[colorz == 2,],2, median,na.rm=TRUE),t='l', col=rgb(0.5,0.5,0.5,1),lwd=2)
points(x = c(42, 56, 70), y = apply(nfatlean[colorz == 3,],2, median,na.rm=TRUE),t='l', col=rgb(0,0,1,0.5),lwd=2)
axis(2, at=c(0,0.2,0.4,0.6), c("0", "20","40", "60"),las=2)
legend("topleft", c("BFMI","Heterozygous", "B6N"), fill = c(rgb(1,0.4,0,0.5),rgb(0.5,0.5,0.5,1),rgb(0,0,1,0.5)), title="UNC5048297 genotype",bg="white",cex=1.2)
axis(1, at= seq(14, 72, 14), seq(14, 72,14) / 7)




cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

lm(as.numeric(growth[,1]) ~ fam + wsize + wlabel + season)

## Litter # + Season + (Family and Litter size per Geno)

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wlabel + season); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(cBfmi <- apply(ngrowth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1] + wsize[colorz==1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(1,0,0,0.5), main = "~ Family(Geno) + litter # + season + Litter size(Geno)")
boxplot(cHetro <- apply(ngrowth[colorz == 2,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 2]+ wsize[colorz==2]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(cB6n <- apply(ngrowth[colorz == 3,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 3]+ wsize[colorz==3]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(cHetro),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(cB6n),d=2), "\n")

lm(as.numeric(ngrowth[colorz == 1,1]) ~ fam[colorz == 1] + wsize[colorz==1])
lm(as.numeric(ngrowth[colorz == 2,1]) ~ fam[colorz == 2]+ wsize[colorz==2])
lm(as.numeric(ngrowth[colorz == 3,1]) ~ fam[colorz == 3]+ wsize[colorz==3])


## Litter # + Season + (Family and Litter size per Geno)

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wlabel + season); return(aa$coefficients["(Intercept)"] + aa$residuals) })

cBfmi <- apply(ngrowth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1] + wsize[colorz == 1]); return(aa$coefficients["(Intercept)"] + aa$residuals) })
rownames(cBfmi) <- rownames(growth[colorz == 1,])
merged <- apply(ngrowth[colorz != 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz != 1] + wsize[colorz != 1]); return(aa$coefficients["(Intercept)"] + aa$residuals) })
rownames(merged) <- rownames(growth[colorz != 1,])

boxplot(cBfmi, col=rgb(1,0,0,0.5), main = "~ Family(Geno) + litter # + season + Litter size(Geno)")
boxplot(merged[colorz[colorz != 1] == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(merged[colorz[colorz != 1] == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(merged[colorz[colorz != 1] == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(merged[colorz[colorz != 1] == 3,]),d=2), "\n")

alldata <- rbind(cBfmi, merged)

result <- matrix(NA, nrow(genotypes[,rownames(alldata)]), ncol(alldata))
for(x in 1:nrow(genotypes[,rownames(alldata)])){
  for(y in 1:ncol(alldata)){
    tryCatch(res <- anova(lm(alldata[,y] ~ as.factor(as.character(genotypes[x, rownames(alldata)]))))[[5]], 
        error = function(e){ res <<- rep(NA, 2) })
    result[x, y] <- res[1]
  }
}


lm(as.numeric(ngrowth[colorz == 1,1]) ~ fam[colorz == 1] + wsize[colorz==1])
lm(as.numeric(ngrowth[colorz != 1,1]) ~ fam[colorz != 1] + wsize[colorz!=1])




### ONLY H ###

genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes  <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                         # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]                                                 # The F2 individuals

fatherGeno <- unlist(genotypes["UNC5048297", as.character(phenotypes[F2, c("Vater")])])
mutterGeno <- unlist(genotypes["UNC5048297", as.character(phenotypes[F2, c("Mutter")])])

F2 <- F2[which(fatherGeno == "H" & mutterGeno == "H")]


colorz <- as.numeric(as.factor(unlist(genotypes["UNC5048297", F2])))
growth <- phenotypes[F2, c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")]

vater <- as.factor(phenotypes[F2, c("Vater")]) ; mutter <- as.factor(phenotypes[F2, c("Mutter")]) ; wsize <- as.numeric(phenotypes[F2, c("WG2")]) ; wlabel <- as.factor(phenotypes[F2, c("W.Label")]) ; season <- as.factor(phenotypes[F2, c("Season")])
fam <- as.factor(paste0(mutter, vater))

## Current model

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ vater + wlabel + season + wsize); return(aa$coefficients["(Intercept)"] + aa$residuals) })
boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Family + Litter size + Litter # + season")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

lm(as.numeric(growth[,1]) ~ fam + wsize + wlabel + season)


## Litter # + Season + (Family and Litter size per Geno)

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wlabel + season); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(cBfmi <- apply(ngrowth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1] + wsize[colorz==1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(1,0,0,0.5), main = "~ Family(Geno) + litter # + season + Litter size(Geno)")
boxplot(cHetro <- apply(ngrowth[colorz == 2,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 2]+ wsize[colorz==2]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(cB6n <- apply(ngrowth[colorz == 3,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 3]+ wsize[colorz==3]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Heterozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(cHetro),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(cB6n),d=2), "\n")

lm(as.numeric(ngrowth[colorz == 1,1]) ~ fam[colorz == 1] + wsize[colorz==1])
lm(as.numeric(ngrowth[colorz == 2,1]) ~ fam[colorz == 2]+ wsize[colorz==2])
lm(as.numeric(ngrowth[colorz == 3,1]) ~ fam[colorz == 3]+ wsize[colorz==3])












lm(as.numeric(growth[,1]) ~ fam + wsize + wlabel + season)


boxplot(apply(ngrowth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz == 1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="red")
boxplot(apply(ngrowth[colorz == 2,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz == 2]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="blue", add=TRUE)
boxplot(apply(ngrowth[colorz == 3,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz == 3]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="blue", add=TRUE)




lm(growth[,"d70"] ~ fam + wsize + wlabel + season + as.factor(colorz))
lm(ngrowth[,"d70"] ~ as.factor(colorz))

op <- par(mfrow=c(1,3))
boxplot(ngrowth[colorz == 1, "d70"], col="red")
boxplot(ngrowth[colorz == 2, "d70"], col="gray")
boxplot(ngrowth[colorz == 3, "d70"], col="blue")



ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ vater + season); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(apply(ngrowth[colorz==1,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz==1] + wlabel[colorz==1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="red")
boxplot(apply(ngrowth[colorz!=1,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz!=1] + wlabel[colorz!=1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="blue", add=TRUE)
boxplot(apply(ngrowth[colorz!=1,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz!=1] + wlabel[colorz!=1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="blue", add=TRUE)




anova(lm(as.numeric(wsize) ~ as.factor(fatherGeno)))

anova(lm(as.numeric(growth[,3]) ~ vater  + wlabel + season + wsize  + colorz:wsize))
anova(lm(as.numeric(wsize) ~ vater + colorz))

lm(as.numeric(ngrowth[colorz==1,3]) ~ wsize[colorz==1])
lm(as.numeric(ngrowth[colorz==2,3]) ~ wsize[colorz==2])
lm(as.numeric(ngrowth[colorz==3,3]) ~ wsize[colorz==3])



anova(lm(as.numeric(growth[colorz==1,3]) ~ vater[colorz==1] + wsize[colorz==1] + wlabel[colorz==1] + season[colorz==1]))
anova(lm(as.numeric(growth[colorz==2,3]) ~ vater[colorz==2] + wsize[colorz==2] + wlabel[colorz==2] + season[colorz==2]))
anova(lm(as.numeric(growth[colorz==3,3]) ~ vater[colorz==3] + wsize[colorz==3] + wlabel[colorz==3] + season[colorz==3]))

### ONLY H

fatherGeno <- unlist(genotypes["UNC5048297", as.character(phenotypes[F2, c("Vater")])])
mutterGeno <- unlist(genotypes["UNC5048297", as.character(phenotypes[F2, c("Mutter")])])

F2 <- F2[which(fatherGeno == mutterGeno & fatherGeno == "H")]

colorz <- as.numeric(as.factor(unlist(genotypes["UNC5048297", F2])))
growth <- phenotypes[F2, c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")]

vater <- as.factor(phenotypes[F2, c("Vater")])
wsize <- as.numeric(phenotypes[F2, c("WG2")])
wlabel <- as.factor(phenotypes[F2, c("W.Label")])
season <- as.factor(phenotypes[F2, c("Season")])

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ vater + wlabel + season + wsize); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(ngrowth[colorz == 1, ], col="red", main = "only H")
boxplot(ngrowth[colorz != 1, ], col="blue", add=TRUE)

lm(as.numeric(growth[,3]) ~ vater + wlabel + season + wsize )

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ vater + wlabel + season); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(apply(ngrowth[colorz==1,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz==1] ); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="red", main = "only H")
boxplot(apply(ngrowth[colorz==2,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz==2] ); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="gray", add=TRUE)
boxplot(apply(ngrowth[colorz==3,],2,function(x){ aa <- lm(as.numeric(x) ~ wsize[colorz==3] ); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col="blue", add=TRUE)

lm(as.numeric(growth[colorz==1,3]) ~ vater + wsize )

lm(as.numeric(ngrowth[colorz==1,3]) ~ wsize[colorz==1] )
lm(as.numeric(ngrowth[colorz==2,3]) ~ wsize[colorz==2] )
lm(as.numeric(ngrowth[colorz==3,3]) ~ vater[colorz==3] + wsize[colorz==3] )

plot(c(1, 9), c(0,max(growth,na.rm=TRUE)), t='n', ylab="Bodyweight")
cnt <- 1
apply(growth, 1, function(x){ points(x, t='l', col = colorz[cnt]); cnt <<- cnt + 1 })


mri <- phenotypes[F2, c("mri42d_fat", "mri42d_lean", "mri56d_fat", "mri56d_lean", "mri70d_fat", "mri70d_lean")]
mri[which(mri[,5] < 0), ]                                               # Which individuals have negative fat weight ?
# 6661341

# Fat mass
plot(c(1, 3), c(0,max(mri[,c(1, 3, 5)], na.rm=TRUE)), t='n', ylab="FAT")
apply(mri[,c(1, 3, 5)],1,function(x){ points(x, t='l', col=2+ sign(cor(1:3, as.numeric(x)))) })
updown <- apply(mri[,c(1, 3, 5)],1,function(x){ sign(cor(1:3, as.numeric(x))) })
#661105 6661107 6661108 6661124 6661156 6661227 6661261 6661308 6661326 6661340 6661341 6661342 6661343 6661386 6661389 6661408 6661448 6661456 6661457 6661460 6661554 6661561 6661758 6661817 
#6662171 6662172 6662173 6662243 6662246 


# Lean mass
plot(c(1, 3), c(0,max(mri[,c(2, 4, 6)], na.rm=TRUE)), t='n', ylab="LEAN")
apply(mri[,c(2, 4, 6)],1,function(x){ points(x, t='l', col=2+ sign(cor(1:3, as.numeric(x)))) })
updown <- apply(mri[,c(2, 4, 6)],1,function(x){ sign(cor(1:3, as.numeric(x))) })
# 6661339 6661341 

