# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/DNA/MegaMuga/")                                                                               # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)                              # Normal A, H, B genotypes

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),phenos]                   # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F2 <- F2[-which(F2=="6661459")]     

topMarker <- "UNC5048297"

inThe1Group <- which(as.character(genotypes[topMarker,F2]) == 1)
outGroup <- which(as.character(genotypes[topMarker,F2]) != 1)

cat("Mean fat levels: AA / H and B: ", mean(as.numeric(phenotypes[F2,"mri56d_fat"])[inThe1Group], na.rm=TRUE), "/", mean(as.numeric(phenotypes[F2,"mri56d_fat"])[outGroup], na.rm=TRUE),"\n")

genotypesIn1 <- genotypes[,F2[inThe1Group]]

#FATres <- lm(as.numeric(phenotypes[F2,"mri42d_fat"]) ~ as.numeric(phenotypes[F2,"WG2"]))$residuals

lowFAT  <- which(as.numeric(phenotypes[F2,"mri42d_fat"])[inThe1Group] < 3.5)
highFAT <- which(as.numeric(phenotypes[F2,"mri42d_fat"])[inThe1Group] > 10)

cat("Group sizes low FAT/high FAT:", length(lowFAT), "/", length(highFAT),"\n")

lowFATtable  <- apply(genotypesIn1[,lowFAT],  1, table)
highFATtable <- apply(genotypesIn1[,highFAT], 1, table)

lowFATnames  <- lapply(lowFATtable, names)
lowFATngeno  <- unlist(lapply(lowFATtable, length))

lowFAT1geno  <- lowFATnames[which(lowFATngeno == 1)]
highFATtable <- highFATtable[names(lowFAT1geno)]

highFATtopgeno       <- lapply(highFATtable, which.max)
highFATtopgenotypes  <- unlist(lapply(highFATtopgeno, names))

LowHighdiff <- names(which(highFATtopgenotypes != unlist(lowFAT1geno)))

plot(rownames(genotypes) %in% LowHighdiff, col=as.numeric(as.factor(map[,"Chr"])), pch=19, cex=0.7)
map[LowHighdiff,]
