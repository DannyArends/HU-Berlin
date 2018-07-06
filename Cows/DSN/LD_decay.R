#LD in Holstein

library(heterozygous)

setwd("~/PopStructCow")

genotypes <- read.table("~/NAS/Cattle/DSN/analysis/QTL/genotypes.txt", sep="\t", check.names=FALSE)
phenotypes <- read.table("~/NAS/Cattle/DSN/analysis/QTL/phenotypes.txt", sep="\t")
topMarkers <- read.table("~/NAS/Cattle/DSN/analysis/QTL/signifcant_marker.txt", sep="\t", check.names=FALSE,header=TRUE)

topMarkers.HF <- topMarkers[topMarkers$breed == "Holstein",]
HFind <- rownames(phenotypes)[which(phenotypes[,"breed"] == "Holstein")]

genotypes <- genotypes[, HFind]
dim(genotypes)
length(HFind)


splitted <- strsplit(rownames(genotypes), "_")
marker.chr <- gsub("Chr", "", unlist(lapply(splitted,"[",1)))
marker.pos <- as.numeric(unlist(lapply(splitted,"[",2)))

Dprimes <- vector("list", nrow(topMarkers.HF))
names(Dprimes) <- paste0("Chr", topMarkers.HF[, "CHR"], "_", topMarkers.HF[, "POS"])
for(x in 1:nrow(topMarkers.HF)){
  chr <- topMarkers.HF[x, "CHR"]
  pos <- as.numeric(topMarkers.HF[x, "POS"])
  top.snp <- paste0("Chr", chr, "_", pos)
  top.numvalues <- rep(NA, ncol(genotypes))
  top.numvalues[genotypes[top.snp,] == -1] <- "A"
  top.numvalues[genotypes[top.snp,] == 0] <- "H"
  top.numvalues[genotypes[top.snp,] == 1] <- "B"
  
  snpInRange <- which(marker.chr == chr & marker.pos >= (pos-5000000) & marker.pos <= (pos+5000000))
  lds <- rep(NA, length(snpInRange))
  names(lds) <- rownames(genotypes)[snpInRange]
  for(snp in snpInRange){
    snp.numvalues <- rep(NA, ncol(genotypes))
    snp.numvalues[genotypes[snp,] == -1] <- "A"
    snp.numvalues[genotypes[snp,] == 0] <- "H"
    snp.numvalues[genotypes[snp,] == 1] <- "B"
    ld <- calculateLD(as.character(top.numvalues), as.character(snp.numvalues))
    lds[rownames(genotypes)[snp]] <- ld$"D'"
  }
  cat("Done", top.snp, "\n")
  Dprimes[[x]] <- lds
}
save(Dprimes, file="~/PopStructCow/linkage.HF.Rdata")

# LD in DSN
library(heterozygous)

setwd("~/PopStructCow")

genotypes <- read.table("~/NAS/Cattle/DSN/analysis/QTL/genotypes.txt", sep="\t", check.names=FALSE)
phenotypes <- read.table("~/NAS/Cattle/DSN/analysis/QTL/phenotypes.txt", sep="\t")
topMarkers <- read.table("~/NAS/Cattle/DSN/analysis/QTL/signifcant_marker.txt", sep="\t", check.names=FALSE,header=TRUE)

topMarkers.DSN <- topMarkers[topMarkers$breed == "DSN",]
DSNind <- rownames(phenotypes)[which(phenotypes[,"breed"] == "DSN")]

genotypes <- genotypes[, DSNind]
dim(genotypes)
length(DSNind)

splitted <- strsplit(rownames(genotypes), "_")
marker.chr <- gsub("Chr", "", unlist(lapply(splitted,"[",1)))
marker.pos <- as.numeric(unlist(lapply(splitted,"[",2)))

Dprimes <- vector("list", nrow(topMarkers.DSN))
names(Dprimes) <- paste0("Chr", topMarkers.DSN[, "CHR"], "_", topMarkers.DSN[, "POS"])
for(x in 1:nrow(topMarkers.DSN)){
  chr <- topMarkers.DSN[x, "CHR"]
  pos <- as.numeric(topMarkers.DSN[x, "POS"])
  top.snp <- paste0("Chr", chr, "_", pos)
  top.numvalues <- rep(NA, ncol(genotypes))
  top.numvalues[genotypes[top.snp,] == -1] <- "A"
  top.numvalues[genotypes[top.snp,] == 0] <- "H"
  top.numvalues[genotypes[top.snp,] == 1] <- "B"
  
  snpInRange <- which(marker.chr == chr & marker.pos >= (pos-5000000) & marker.pos <= (pos+5000000))
  lds <- rep(NA, length(snpInRange))
  names(lds) <- rownames(genotypes)[snpInRange]
  for(snp in snpInRange){
    snp.numvalues <- rep(NA, ncol(genotypes))
    snp.numvalues[genotypes[snp,] == -1] <- "A"
    snp.numvalues[genotypes[snp,] == 0] <- "H"
    snp.numvalues[genotypes[snp,] == 1] <- "B"
    ld <- calculateLD(as.character(top.numvalues), as.character(snp.numvalues))
    lds[rownames(genotypes)[snp]] <- ld$"D'"
  }
  cat("Done", top.snp, "\n")
  Dprimes[[x]] <- lds
}
save(Dprimes, file="~/PopStructCow/linkage.DSN.Rdata")