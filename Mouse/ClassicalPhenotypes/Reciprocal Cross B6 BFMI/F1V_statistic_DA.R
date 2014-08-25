# F1V_statistic_DA.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Stefan Kärst, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")
f <- read.table(file="20140718_F1V.txt", sep="\t", na.strings = "NA", stringsAsFactors=F, header=TRUE )

mb6 <- subset(f, f$dirocross == "matB6" & f$futter == "FF"); cat("Number of animals (m+f, maternalB6, Fat Food):", nrow(mb6), "\n")
mbf <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF"); cat("Number of animals (m+f, maternalBFMI, Fat Food):", nrow(mbf), "\n")
mb6m <- subset(f, f$dirocross == "matB6" & f$futter == "FF" & f$sex== "m"); cat("Number of animals (male, maternalB6, Fat Food):", nrow(mb6m), "\n")
mbfm <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF" & f$sex== "m"); cat("Number of animals (male, maternalBFMI, Fat Food):", nrow(mbfm), "\n")
mb6f <- subset(f, f$dirocross == "matB6" & f$futter == "FF" & f$sex== "f"); cat("Number of animals (female, maternalB6, Fat Food):", nrow(mb6f), "\n")
mbff <- subset(f, f$dirocross == "matBFMI" & f$futter == "FF" & f$sex== "f"); cat("Number of animals (female, maternalBFMI, Fat Food):", nrow(mbff), "\n")


plotStefan <- function(subsetA, subsetB, what="fat2lean", yMin = 0, yMax= 0.31){
  y <- apply(subsetA[,paste("d", x, what, sep="")], 2, mean)
  ySD <- apply(subsetA[,paste("d", x, what, sep="")], 2, sd)

  plot(x = c(28,70), y = c(yMin, yMax), type='n', xlab="Age [days]", ylab="Fat Mass/Lean Mass Ratio", main=paste("Fat/Lean Mass Ratio\n", "F1 from reciprocal BFMI x B6, under HFD"))
  grid(col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
  points(x, y, cex=2, type='l', col=B6)
  plotCI(x, y, uiw=ySD, add=T, gap=0, col=B6, pch=19, cex=2)

  # Fat Food, matBFMI
  y <- apply(subsetB[,paste("d", x, what, sep="")], 2, mean)
  ySD <- apply(subsetB[,paste("d", x, what, sep="")], 2, sd)
  plotCI(x, y, uiw=ySD, add=T, col=BFMI, gap=0, pch=16, cex=2)
  points(x, y, cex=2, type='l', col=BFMI)
  legend("bottomright", pch=c(16,16), col=c(B6,BFMI), c("m+f, matB6", "m+f, matBFMI"))
}

#### Plot 1: Fat/Lean
plotStefan(mb6, mbf)

plotStefan(mb6, mbf, "", 0, 40)

plotStefan(mb6m, mbfm)

plotStefan(mb6m, mbfm, "", 0, 40)

plotStefan(mb6f, mbff)

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

