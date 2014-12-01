# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

getSeason <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES),".", fixed=TRUE),"[",2)))
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5] <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8] <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11] <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2] <- "Winter"
  return(ret)
}

setwd("E:/Mouse/DNA/MegaMuga/")                                                                               # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE)                              # Normal A, H, B genotypes

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("Vater", "W.dat", "W.Label", "Eltern_ID", "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),phenos]                   # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F2 <- F2[-which(F2=="6661459")]     

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix

onegenotype <- which(lapply(apply(genotypes[,F2], 1, table), length) == 1)                                    # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- unique(genotypes[-onegenotype, F2])                                                            # Only take the F2 individuals, and the unique markers
cat("Left with", nrow(genotypes), "markers\n")                                                                # == Left with 11677 markers

topMarker <- "UNC5048297"
inThe1Group <- which(as.character(genotypes[topMarker,F2]) == 1)
outGroup <- which(as.character(genotypes[topMarker,F2]) != 1)

cat("Mean fat levels: AA / H and B: ", mean(as.numeric(phenotypes[F2,"mri56d_fat"])[inThe1Group], na.rm=TRUE), "/", mean(as.numeric(phenotypes[F2,"mri56d_fat"])[outGroup], na.rm=TRUE),"\n")

genotypesIn1 <- genotypes[,F2[inThe1Group]]

#FATres <- lm(as.numeric(phenotypes[F2,"mri42d_fat"]) ~ as.numeric(phenotypes[F2,"WG2"]))$residuals

subfamily      <- as.factor(phenotypes[F2, "Vater"])                                                                     # Fixed effect: Subfamily structure (factor)
littersize     <- as.numeric(phenotypes[F2, "WG2"])                                                                      # Fixed effect: Size of the litter  (linear effect)
litternumber   <- as.factor(phenotypes[F2, "W.Label"])                                                                   # Fixed effect: Number of litter    (factor)
season         <- as.factor(phenotypes[F2, "Season"])                                                                    # Fixed effect: Season when born    (factor)

myresiduals <- lm(phenotypes[F2,"mri42d_fat"] ~ subfamily + littersize + litternumber + season)$residuals

newFatPhenotype <- phenotypes[F2,"mri42d_fat"]            # get the original phenotype
names(newFatPhenotype) <- 1:length(newFatPhenotype)       # names
newFatPhenotype[names(myresiduals)] <- myresiduals        # Fill in the residuals

newPhenotype <- phenotypes[F2,"mri42d_fat"]            # get the original phenotype
names(newPhenotype) <- 1:length(newPhenotype)       # names
newPhenotype[names(aa$fitted.values )] <- aa$fitted.values        # Fill in the residuals

lowFAT  <- which(as.numeric(newFatPhenotype)[inThe1Group] < 0)
highFAT <- which(as.numeric(newFatPhenotype)[inThe1Group] > 0)

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
