
setwd("E:/Mouse/DNA/MegaMuga/")

# Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                              # Normal A, H, B genotypes
genotypesPh <- read.table("Analysis/genotypesPhasedBeagle.txt", sep="\t", check.names=FALSE, colClasses="character", na.strings=c("","AA","CC","TT","GG"))        # Phased by Beagle (only heterozygous)
genotypesGP <- read.table("Analysis/genotypesPhasedGP.txt",     sep="\t", check.names=FALSE, colClasses="character", na.strings=c(""))                            # Phased towards the grandparents

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                                                                    # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("Vater", "W.dat", "W.Label", "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
            
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)), phenos]               # Use only the phenotypes for which we have genotypes

setwd("E:/Mouse/DNA/MegaMuga/")

toRem <- which(rownames(phenotypes) %in% c("6661459","6661341", "6661339"))                                # Remove 3 strange individuals from the phenotypes (1x no genotypes, 2x weird mass measurements)

write.table(phenotypes[-toRem,], file="Phenotypes/MatchedPhenotypes.txt", sep="\t")