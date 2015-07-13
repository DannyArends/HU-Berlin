#
# New QTL mapping of BFMI, using a slightly different population structure correction
#

source("D:/Github/HU-Berlin/Mouse/Muga/dateToSeason.R")
setwd("E:/Mouse/DNA/MegaMuga/")                                                                                           # Read in the data from the MegaMuga

map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")      # Normal A, H, B genotypes

missingPerInd <- apply(genotypes, 2, function(x){ sum(is.na(x)) / length(x) * 100 })
wrongIND <- c("6661965", "6662155", "6662156", names(which(missingPerInd==100)))

phenotypes  <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                                 # Add the season column to the matrix
phenotypes <- phenotypes[-which(rownames(phenotypes) %in% wrongIND),]                                                     # Remove the faulty individuals from the phenotypes
genotypes <- genotypes[,-which(colnames(genotypes) %in% wrongIND)]                                                        # Remove the faulty individuals from the genotypes
            
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                             # The F2 individuals
F1 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 27)]                                                             # The F2 individuals
            
genoDad <- genotypes[, as.character(phenotypes[F2,"Vater"])] ; colnames(genoDad) <- F2            
onegenotype <- which(lapply(apply(genoDad, 1, table), length) < 2)                                                        # Markers with only one genotype cannot be used in QTL mapping
genoDad <- genoDad[-onegenotype,]           
            
genoMom <- genotypes[, as.character(phenotypes[F2,"Mutter"])] ; colnames(genoMom) <- F2           
onegenotype <- which(lapply(apply(genoMom, 1, table), length) < 2)                                                        # Markers with only one genotype cannot be used in QTL mapping
genoMom <- genoMom[-onegenotype,]

idx <- which(phenotypes[F2,"W.Label"] == "A")

resWGdad <- apply(genoDad, 1, function(x){
                                          tryCatch(
                                              res <- anova(lm(phenotypes[F2,"WG"][idx] ~ x[idx]))[[5]], 
                                              error = function(e){ res <<- rep(NA, 2) })
                                            return(res) })

