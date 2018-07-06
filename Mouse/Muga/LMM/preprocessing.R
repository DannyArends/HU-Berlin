# Preprocessing of the MegaMuga data
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

setwd("D:/Edrive/Mouse/DNA/MegaMuga/") # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE, colClasses="character") # Normal A, H, B genotypes

parentgeno <- genotypes[,c("BFMI860-12 (V2)", "B6N")]

setwd("D:/Edrive/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI") # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("Vater", "W.dat", "W.Label", "Eltern_ID", "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),phenos] # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)] # The F2 individuals
F2 <- F2[-which(F2 %in% c("6661339", "6661341","6661459", "6661965", "6662155", "6662156"))]

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"])) # Add the season column to the matrix
m <- "JAX00191930"

genotypes.all <- genotypes
genotypes <- genotypes[,F2]
for(m in 1:nrow(genotypes)){
  mtable <- table(as.character(genotypes[m,]))
  small <- which(mtable < 30)
  if(length(small) > 0) genotypes[m, which(genotypes[m, ] %in% names(small))] <- NA
  if(m %% 1000 == 0)cat(m, "/", nrow(genotypes), "\n")
}

onegenotype <- which(lapply(apply(genotypes, 1, table), length) <= 1) # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- genotypes[-onegenotype,] # Only take segregating markers
cat("Left with", nrow(genotypes), "markers\n") # == Left with 10949 markers

founder.geno <- genotypes.all[rownames(genotypes), c("BFMI860-12 (V2)", "B6N")]
genotypes.all <- genotypes.all[rownames(genotypes), F2]
phenotypes <- phenotypes[F2,]

founder.issues <- unique(c(which(founder.geno[,1] == "H"), which(founder.geno[,2] == "H"), # Founders cannot be Heterozygous
                           which(apply(founder.geno,1,function(x){any(is.na(x))})),  # Founders cannot be NA
                           which(apply(founder.geno,1,function(x){return(x[1] == x[2])})))) # Founder cannot be equal

founder.geno <- founder.geno[-founder.issues,]
genotypes.all <- genotypes.all[rownames(founder.geno), ]
genotypes <- genotypes[rownames(founder.geno), ]
map <- map[rownames(genotypes),]
cat("Left with", nrow(genotypes), "markers\n") # == Left with 10949 markers

m <- 1
genotypes.recoded <- genotypes
for(mname in rownames(genotypes)){
  genotypes.recoded[mname, which(genotypes.recoded[mname,] == founder.geno[mname, 2])] <- "N"
  genotypes.recoded[mname, which(genotypes.recoded[mname,] == founder.geno[mname, 1])] <- "B"
  m <- m + 1
  if(m %% 1000 == 0)cat(m, "/", nrow(genotypes), "\n")
}


map <- cbind(map, sumPos = NA)
chrs <- 1:20
names(chrs) <- c(1:19, "X")

chrs.starts <- c(0)
chrs.lengths <- c()
chrs.summed <- 0
chr.gap <- 25000000

for(chr in names(chrs)){
  onChr <- which(map[,"Chr"] == chr)
  chr.length <- max(as.numeric(map[onChr, "Mb_NCBI38"]))
  map[onChr,"sumPos"] <- as.numeric(map[onChr,"Mb_NCBI38"]) + chrs.summed
  chrs.summed <- chrs.summed + chr.length + chr.gap
  chrs.lengths <- c(chrs.lengths, chr.length + chr.gap)
  chrs.starts <- c(chrs.starts, chrs.summed)
}
chrs.lengths <- c(chrs.lengths, NA)

setwd("D:/Edrive/Mouse/DNA/MegaMuga/inputF2")
write.table(genotypes, file="cleaned_genotypes_F2.txt", sep='\t', quote=FALSE)
write.table(genotypes.recoded, file="cleaned_recoded_genotypes_F2.txt", sep='\t', quote=FALSE)
write.table(genotypes.all, file="all_genotypes_F2.txt", sep='\t', quote=FALSE)
write.table(founder.geno, file="cleaned_genotypes_founders.txt", sep='\t', quote=FALSE)
write.table(map, file="cleaned_map_25MbGap.txt", sep='\t', quote=FALSE)
write.table(phenotypes, file="cleaned_phenotypes_F2.txt", sep='\t', quote=FALSE)
