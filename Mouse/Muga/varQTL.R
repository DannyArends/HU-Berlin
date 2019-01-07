 library(plotrix)

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")                                                                               # Read in the data from the Mega Muga
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE, colClasses="character")      # Normal A, H, B genotypes
phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)

F2 <- c(rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)])                                                # The F2 individuals

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix

gt <- factor(unlist(genotypes["UNC5048297", F2]), levels = c("A", "H", "B"))

indAA <- names(gt)[which(gt == "A")]
indH <- names(gt)[which(gt == "H")]
indBB <- names(gt)[which(gt == "B")]


pheAA <- phenotypes[indAA, paste0("d", seq(21, 70,7))]
pheH <- phenotypes[indH, paste0("d", seq(21, 70,7))]
pheBB <- phenotypes[indBB, paste0("d", seq(21, 70,7))]

statsA <- apply(pheAA, 2, function(x){ return(c(mean(x), sd(x)))})
statsH <- apply(pheH, 2, function(x){ return(c(mean(x), sd(x)))})
statsB <- apply(pheBB, 2, function(x){ return(c(mean(x), sd(x)))})

statsA <- rbind(statsA, statsA[2,] / statsA[1,])
statsH <- rbind(statsH, statsH[2,] / statsH[1,])
statsB <- rbind(statsB, statsB[2,] / statsB[1,])

boxplot(cbind(statsA[3,], statsH[3,], statsB[3,]))


t.test(statsA[2,],statsB[2,])