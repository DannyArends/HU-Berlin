# F1V_statistic_DA.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Stefan Kärst, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

library(plotrix)

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")
f <- read.table(file="20140718_F1V.txt", sep="\t", na.strings = "NA", stringsAsFactors=F, header=TRUE )

mb6 <- subset(f, f$dirocross == "matB6" & f$futter == "FF"); cat("Number of animals (m+f, maternalB6, Fat Food):", nrow(mb6), "\n")
mbf <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF"); cat("Number of animals (m+f, maternalBFMI, Fat Food):", nrow(mbf), "\n")
mb6_ZF <- subset(f, f$dirocross == "matB6" & f$futter == "ZF"); cat("Number of animals (m+f, maternalB6, Fat Food):", nrow(mb6), "\n")
mbf_ZF <- subset(f, f$dirocross == "matBFMI" & f$futter == "ZF"); cat("Number of animals (m+f, maternalBFMI, Fat Food):", nrow(mbf), "\n")

mb6m <- subset(f, f$dirocross == "matB6" & f$futter == "FF" & f$sex== "m"); cat("Number of animals (male, maternalB6, Fat Food):", nrow(mb6m), "\n")
mbfm <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF" & f$sex== "m"); cat("Number of animals (male, maternalBFMI, Fat Food):", nrow(mbfm), "\n")
mb6m_ZF <- subset(f, f$dirocross == "matB6" & f$futter == "ZF" & f$sex== "m"); cat("Number of animals (male, maternalB6, Fat Food):", nrow(mb6m), "\n")
mbfm_ZF <- subset(f, f$dirocross == "matBFMI" & f$futter == "ZF" & f$sex== "m"); cat("Number of animals (male, maternalBFMI, Fat Food):", nrow(mbfm), "\n")


mb6f <- subset(f, f$dirocross == "matB6" & f$futter == "FF" & f$sex== "f"); cat("Number of animals (female, maternalB6, Fat Food):", nrow(mb6f), "\n")
mbff <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF" & f$sex== "f"); cat("Number of animals (female, maternalBFMI, Fat Food):", nrow(mbff), "\n")

B6 <- "gray" ; BFMI <- "orange"

plotStefan <- function(subsetA, subsetB, what="fat2lean", yMin = 0, yMax= 0.31, xRange  = seq(28,70,7), legendtext =  c("m+f, matB6", "m+f, matBFMI"), ...){
  y <- apply(subsetA[,paste("d", xRange, what, sep="")], 2, mean)
  ySD <- apply(subsetA[,paste("d", xRange, what, sep="")], 2, sd)

  plot(x = c(min(xRange),70), y = c(yMin, yMax), type='n', ...)
  grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
  points(xRange, y, cex=2, type='l', col=B6)
  plotCI(xRange, y, uiw=ySD, add=T, gap=0, col=B6, pch=19, cex=2)

  # Fat Food, matBFMI
  y <- apply(subsetB[,paste("d", xRange, what, sep="")], 2, mean)
  ySD <- apply(subsetB[,paste("d", xRange, what, sep="")], 2, sd)
  plotCI(xRange, y, uiw=ySD, add=T, col=BFMI, gap=0, pch=16, cex=2)
  points(xRange, y, cex=2, type='l', col=BFMI)
  legend("bottomright", pch=c(16,16), col=c(B6,BFMI), legendtext)
}

#### Plot 1: Fat/Lean (HFD and SBD)
png("Analysis/plot1.png", width=1200, height=600)
  op <- par(mfrow = c(1, 2))
  plotStefan(mb6, mbf, xlab="Age [days]", ylab="Fat Mass/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under HFD"))
  plotStefan(mb6_ZF, mbf_ZF, xlab="Age [days]", ylab="Fat Mass/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under SBD"))
dev.off()

#### Plot 2: Bodyweight (HFD and SBD)
png("Analysis/plot2.png", width=1200, height=600)
  op <- par(mfrow = c(1, 2))
  plotStefan(mb6, mbf, what="", 0, 40, xRange = seq(21,70,7), xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under HFD"))
  plotStefan(mb6_ZF, mbf_ZF, what="", 0, 40, xRange = seq(21,70,7), xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under SBD"))
dev.off()

#### Plot 3: Bodyweight (HFD and SBD)
png("Analysis/plot3.png", width=1200, height=600)
  op <- par(mfrow = c(1, 2))
  plotStefan(mb6m, mbfm, what="", 0, 40, legendtext = c("male, matB6", "male, matBFMI"), xRange = seq(21,70,7), xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under HFD"))
  plotStefan(mb6m_ZF, mbfm_ZF, what="", 0, 40, legendtext = c("male, matB6", "male, matBFMI"), xRange = seq(21,70,7), xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under SBD"))
dev.off()

#### Plot 4: Fat/Lean
plotStefan(mb6m, mbfm, legendtext = c("male, matB6", "male, matBFMI"), xlab="Age [days]", ylab="Fat Mass/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under HFD"))

#### Plot 5: Bodyweight *male only*
plotStefan(mb6m, mbfm, "", 0, 40, legendtext = c("male, matB6", "male, matBFMI"), xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under HFD"))

#### Plot 6: Bodyweight *female only*
plotStefan(mb6f, mbff, "", 0, 40, legendtext = c("female, matB6", "female, matBFMI"), xlab="Age [days]", ylab="Bodyweight", main=paste("Bodyweight\n", "F1 from reciprocal BFMI x B6, under HFD"))

#### Plot 7: Fat/Lean *female only*
plotStefan(mb6f, mbff, legendtext = c("female, matB6", "female, matBFMI"), xlab="Age [days]", ylab="Fat Mass/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under HFD"))

paste("d", c(21, 28, 35, 42, 49, 56, 63, 70))

newmatrix <- f[,c("ID", "vater", "mutter", "dirocross", "WG", "sex", "futter")]
newmatrix <- rbind(newmatrix,newmatrix,newmatrix,newmatrix,newmatrix,newmatrix,newmatrix,newmatrix)
newmatrix <- cbind(newmatrix, day =  unlist(lapply(c(21, 28, 35, 42, 49, 56, 63, 70),rep,60)))
newmatrix <- cbind(newmatrix, bodyweight =  unlist(f[,paste("d", c(21, 28, 35, 42, 49, 56, 63, 70),sep="")]))
newmatrix <- cbind(newmatrix, mrifat    =  c(rep(NA, 60), unlist(f[,paste("d", c(28, 35, 42, 49, 56, 63, 70),"mrifat",sep="")])))
newmatrix <- cbind(newmatrix, mrilean   =  c(rep(NA, 60), unlist(f[,paste("d", c(28, 35, 42, 49, 56, 63, 70),"mrilean",sep="")])))
newmatrix <- cbind(newmatrix, fat2lean  = newmatrix$mrifat / newmatrix$mrilean)
newmatrix <- cbind(newmatrix, mrimass   = newmatrix$mrifat + newmatrix$mrilean)
newmatrix <- cbind(newmatrix, bmi       = newmatrix$mrifat / newmatrix$mrimass)

traits <- c("glukosemg", "glucosemmol", "GF1", "GF2", "totalGF", "RF1", "RF2", "totalRF", "IF", "muskel", "leber", "BAT", "LD")

aaa <- matrix(NA, nrow(newmatrix), length(traits))
colnames(aaa) <- traits
newmatrix <- cbind(newmatrix, aaa)

for(trait in traits){ newmatrix[newmatrix$day==70, trait] <- f[, trait]; }

# Now we have the matrix we can do our testing on:

traits <- c("bodyweight","mrifat", "mrilean", "fat2lean", "mrimass", "bmi")

models <- vector("list", length(traits))
names(models) <- traits
for(trait in traits){
  models[[trait]] <- anova(lm(newmatrix[,trait] ~ as.factor(newmatrix$day) + as.factor(newmatrix$WG) + newmatrix$sex + newmatrix$futter + newmatrix$dirocross))
}

traits <- c("glukosemg", "glucosemmol", "GF1", "GF2", "totalGF", "RF1", "RF2", "totalRF", "IF", "muskel", "leber", "BAT", "LD")
modelsSmall <- vector("list", length(traits))
names(modelsSmall) <- traits
for(trait in traits){
  modelsSmall[[trait]] <- anova(lm(newmatrix[,trait] ~ as.factor(newmatrix$WG) + newmatrix$sex + newmatrix$futter + newmatrix$dirocross))
}

