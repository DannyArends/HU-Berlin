# Muga Phenotypes
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/DNA/MegaMuga/")                                                             # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                              # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),]       # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]
P  <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 26)]

phenos <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG2", "Farbe","sex")

sum(length(F2), length(F1), length(P))

plot(apply(genotypes[which(map[,"Chr"] == "X"), F2], 2, function(x){length(which(x == "H"))}))        # Plot the number of heterozygous loci on the X chromosome
which(apply(genotypes[which(map[,"Chr"] == "X"), F2], 2, function(x){length(which(x == "H"))}) > 20)  # Which animals are more heterozygous then expected
# 6661524 6661965 6662141 6662142 6662143 6662144

map[is.na(map[,"cM"]), "cM"] <- 1                                                                     # make sure there is always a cM position

write.table(rbind(cbind(Chr = "", cM = "", t(phenotypes[F2,phenos])), cbind(map[,c("Chr","cM")],genotypes[,F2])), "cross.csvr",sep=",", quote=FALSE, col.names=FALSE)
library(qtl)
cross <- read.cross("csvr", file="cross.csvr")

Xchr <- genotypes[rownames(map)[which(map$Chr == "X")], F2]                                 # Get the X-chromosomes for the F2 individuals
heterozygous <- apply(Xchr, 2, function(x){ sum(x == "H",na.rm = TRUE) } ) / nrow(Xchr)     # Calculate the amount of heterozygous markers on the X-chromosomes
plot(heterozygous)                                                                          # Plot the amount of heterozygous per individual 

# TODO: Scan the phenotypes using the correct models
# TODO: Add covariates, and interactions
res42 <- scanone(cross, pheno.col="mri42d_fat", addcovar="WG2", model="np")
res56 <- scanone(cross, pheno.col="mri56d_fat", model="np")
res70 <- scanone(cross, pheno.col="mri70d_fat", model="np")
plot(res56)

# JAX00107082 3 47884149-50986795
# UNC5246978 3 50986795