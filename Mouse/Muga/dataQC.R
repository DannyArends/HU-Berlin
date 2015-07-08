# data Quality control
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014
library(lme4)

source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("33310233", "6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes  <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                         # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                            # Remove the faulty individuals from the genotypes

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]                                                 # The F2 individuals

colorz <- as.numeric(as.factor(unlist(genotypes["UNC5048297", F2])))
growth <- phenotypes[F2, c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")]

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
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(growth[colorz == 1,]) - colMeans(growth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(growth[colorz == 1,]) - colMeans(growth[colorz == 3,]),d=2), "\n")

## Family (Fixed)

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ fam ); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Family")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Family (Random)

ngrowth <- apply(growth, 2, function(x){ aa <- lmer(as.numeric(x) ~ (1|fam) ); return(residuals(aa)) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Family")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Litter size

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ wsize ); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Litter size")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Family + Litter size

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ fam + wsize ); return(aa$coefficients["(Intercept)"] + aa$residuals) })

boxplot(ngrowth[colorz == 1,], col=rgb(1,0,0,0.5), main = "~ Family + Litter size")
boxplot(ngrowth[colorz == 2,], col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(ngrowth[colorz == 3,], col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 2,]),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(ngrowth[colorz == 1,]) - colMeans(ngrowth[colorz == 3,]),d=2), "\n")

## Family per Geno

boxplot(cBfmi <- apply(growth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(1,0,0,0.5), main = "~ Family(Geno)")
boxplot(cHetro <- apply(growth[colorz == 2,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 2]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(cB6n <- apply(growth[colorz == 3,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 3]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(cHetro),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(cB6n),d=2), "\n")

## Family and Litter size per Geno

boxplot(cBfmi <- apply(growth[colorz == 1,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 1] + wsize[colorz==1]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(1,0,0,0.5), main = "~ Family(Geno) + Litter size(Geno)")
boxplot(cHetro <- apply(growth[colorz == 2,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 2]+ wsize[colorz==2]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0.5,0.5,0.5,0.5), add=TRUE)
boxplot(cB6n <- apply(growth[colorz == 3,],2,function(x){ aa <- lm(as.numeric(x) ~ fam[colorz == 3]+ wsize[colorz==3]); return(aa$coefficients["(Intercept)"] + aa$residuals) }), col=rgb(0,0,1,0.5), add=TRUE)
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(cHetro),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(cB6n),d=2), "\n")

## Current model

ngrowth <- apply(growth, 2, function(x){ aa <- lm(as.numeric(x) ~ fam + wsize + wlabel + season); return(aa$coefficients["(Intercept)"] + aa$residuals) })

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
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

cat("BFMI - Het:", round(colMeans(cBfmi) - colMeans(cHetro),d=2), "\n")
cat("BFMI - B6n:", round(colMeans(cBfmi) - colMeans(cB6n),d=2), "\n")

lm(as.numeric(ngrowth[colorz == 1,1]) ~ fam[colorz == 1] + wsize[colorz==1])
lm(as.numeric(ngrowth[colorz == 2,1]) ~ fam[colorz == 2]+ wsize[colorz==2])
lm(as.numeric(ngrowth[colorz == 3,1]) ~ fam[colorz == 3]+ wsize[colorz==3])


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
legend("topleft", c("BFMI locus","Hetrozygous", "B6n"), col = c(rgb(1,0,0,0.5),rgb(0.5,0.5,0.5,0.5),rgb(0,0,1,0.5)), lwd=2)

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

