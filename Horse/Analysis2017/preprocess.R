opposite <- function(x){
  ret <- NULL
  for(e in x){
    if(is.na(e)){
      ret <- c(ret, "")
    }else{
      if(e == "A") ret <- c(ret, "T")
      if(e == "C") ret <- c(ret, "G")
      if(e == "G") ret <- c(ret, "C")
      if(e == "T") ret <- c(ret, "A")
    }
  }
  return(ret)
}

# Read the raw data from the excel file which was saved as a TXT file
setwd("D:/Edrive/Horse/DNA/")
arabian  <- read.csv("Equine60k/input/arabianhorses.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=1)

# Load marker information (Illumina)
markerinfo <- read.csv("60Karray/GGP_Equine.csv", skip=7, header=TRUE, colClasses="character")
markerinfo[,"Name"] <- gsub("-", "_", as.character(markerinfo[,"Name"]))
cat("Loaded:", nrow(markerinfo), "markers\n")

# Start the analysis
setwd("D:/Edrive/Horse/DNA/Analysis2017")
genotypes             <- arabian[29:nrow(arabian), ]                  # Row 30 till the end contains genotype data
rownames(genotypes)   <- gsub("-", "_", arabian[29:nrow(arabian), 1])
cat(nrow(genotypes), "SNP calls in the original file\n")
genotypes <- genotypes[which(genotypes[,"GenTrain.Score"] >= 0.6), ]
cat(nrow(genotypes), "SNP calls in the after filtering for GenTrain >= 0.6\n")

map <- genotypes[, 1:3]                  # Row 30 till the end contains genotype map
rownames(map) <- gsub("-", "_", map[, 1])
map <- map[, -1]
write.table(map, "genotypes_map.txt", sep="\t", quote=FALSE, na = "",row.names=TRUE)

genotypes <- genotypes[,c(grep("Top.Alleles", colnames(genotypes)),grep("TOP", colnames(genotypes)))]
colnames(genotypes) <- gsub(".TOP", "", colnames(genotypes))
colnames(genotypes) <- gsub(".Top.Alleles", "", colnames(genotypes))

arabGeno <- matrix(NA, nrow(genotypes), ncol(genotypes), dimnames = list(rownames(genotypes), colnames(genotypes)))
for(x in 1:nrow(arabGeno)){
  mname <- rownames(arabGeno)[x]
  genoIndex <- which(rownames(genotypes) == mname)
  cat(mname, genoIndex, "\n")
  if(length(genoIndex) == 1){
    genos <- unlist(genotypes[genoIndex, ])
    infoIndex <- which(markerinfo[,"Name"] == mname)
    if(length(infoIndex) == 1 && !is.na(markerinfo[infoIndex, "IlmnStrand"])){
      if(markerinfo[infoIndex, "IlmnStrand"] == "BOT" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        arabGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else if(markerinfo[infoIndex, "IlmnStrand"] == "TOP" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        arabGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else{
        arabGeno[mname, ] <- genos
      }
    }else{
      arabGeno[mname, ] <- genos
    }
  }
}
write.table(arabGeno, "genotypes_snp.txt", sep="\t", quote=FALSE, na = "")

phenotypes <- arabian[1:28, ]
phenotypes <- phenotypes[,grep("GType", colnames(phenotypes))]
colnames(phenotypes) <- gsub(".GType", "", colnames(phenotypes))

write.table(phenotypes, "phenotypes.txt", sep="\t", quote=FALSE, na = "")

### Read in the data, order the genotypes by chromosome and location
setwd("D:/Edrive/Horse/DNA/Analysis2017")
genotypes <- read.table("genotypes_snp.txt", sep="\t", na.strings=c("","NA"), colClasses="character", row.names=1, header=TRUE)
map <- read.table("genotypes_map.txt", sep="\t", na.strings=c("","NA"), colClasses=c("character","character","numeric"), row.names=1, header=TRUE)
map <- map[order(map$Position),]

chromosomeOrder <- c(as.character(rep(1:31)), "X", "Y")

orderedMap <- NULL
for(chr in chromosomeOrder){
  orderedMap <- rbind(orderedMap, map[which(map$Chromosome == chr),])
}

genotypes <- genotypes[rownames(orderedMap),]

write.table(orderedMap, "genotypes_map_ordered.txt", sep="\t", quote=FALSE)
write.table(genotypes, "genotypes_snp_ordered.txt", sep="\t", quote=FALSE)

### Re-Read the chromosome ordered data
setwd("D:/Edrive/Horse/DNA/Analysis2017")
phenotypes <- read.table("phenotypes.txt", sep="\t", na.strings="", colClasses="character", row.names=1, header=TRUE)
genotypes <- read.table("genotypes_snp_ordered.txt", sep="\t", na.strings=c("","NA"), colClasses="character", row.names=1, header=TRUE)
map <- read.table("genotypes_map_ordered.txt", sep="\t", na.strings=c("","NA"), colClasses=c("character","character","numeric"), row.names=1, header=TRUE)
all(rownames(genotypes) == rownames(map))

### Filtering down genotypes of individuals older then 4
olderThen4 <- colnames(phenotypes)[which(as.numeric(phenotypes["Age",]) > 4)]       # Only individuals older then 4
genotypes <- genotypes[,olderThen4]
phenotypes <- phenotypes[,olderThen4]

### Filtering non segregating markers
nosegg  <- which(unlist(lapply(apply(genotypes, 1, table),length)) < 2)             # Non informative, since we only have 1 genotype
genotypes <- genotypes[-nosegg,]
map <- map[-nosegg,]

### Filtering markers with too much missing data
toomuchmissing <- which(unlist(apply(genotypes, 1, function(x){ return(sum(is.na(x))) })) / ncol(genotypes) > 0.05)
genotypes <- genotypes[-toomuchmissing,]
map <- map[-toomuchmissing,]

### Filter based on allele frequency
allelefreq <- apply(genotypes, 1 , function(x){
  tbl <- table(unlist(lapply(x, strsplit, "")))
  min(tbl / sum(tbl))
})
toolowallelefreq <- which(allelefreq <= 0.05)

genotypes <- genotypes[-toolowallelefreq,]
map <- map[-toolowallelefreq,]

## Remove markers which have only a couple of individuals in a certain genotype
tables <- apply(genotypes, 1, table)
minN <- lapply(tables, min)
lowN <- which(unlist(minN) <= 3)

genotypes <- genotypes[-lowN,]
map <- map[-lowN,]

### Not in HWE
library(heterozygous)
hwe <- HWE(genotypes)
hweadj <- p.adjust(hwe, "BH")
outofhwe <- which(hweadj < 0.1)

genotypes <- genotypes[-outofhwe,]
map <- map[-outofhwe,]

### Prune using pruneLD.R

write.table(map.pruned, "genotypes_map_filtered.txt", sep="\t", quote=FALSE)
write.table(genotypes.pruned, "genotypes_snp_filtered.txt", sep="\t", quote=FALSE)


### Adjust the phenotype data based on sex
phenotypes <- phenotypes[,olderThen4]
morphological <- c("Strain", "WH", "CW", "CH", "NG", "TG", "ChG", "ChD", "ChW", "BLL", "BL", "FCL", "HCL")
sex <- as.factor(as.character(phenotypes["Sex", ]))

phenotypes.corrected <- matrix(NA, length(morphological), ncol(phenotypes), dimnames = list(morphological, colnames(phenotypes)))

for(trait in morphological) {
  model <- lm(as.numeric(phenotypes[trait, ]) ~ sex)
  newpheno <- model$residuals + model$coefficients["(Intercept)"]
  phenotypes.corrected[trait, ] <- round(newpheno, d = 2)
}

write.table(phenotypes.corrected, "phenotypes_corrected_filtered.txt", sep="\t", quote=FALSE)

#CC TC TT 
#10 18  9
