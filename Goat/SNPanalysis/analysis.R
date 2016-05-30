# Analysis of the SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016
library(ape)

setwd("E:/Goat/DNA/SihamAnalysis")

snpdata    <- read.table("filtered_snps.txt", sep="\t", check.names=FALSE, colClasses="character")
snpinfo    <- read.table("snpinfo.txt", sep="\t", na.strings=c("", "NA", "N.D."))
samples    <- read.table("sampleinfo.txt", sep="\t")
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`
samples    <- cbind(samples, locations[rownames(samples),])
samples    <- cbind(samples, locationShort = as.character(unlist(lapply(strsplit(as.character(samples[,"Location"]), "_"), "[",1))))

snpAlleles <- lapply(strsplit(as.character(snpinfo[,"allele"]), ""), "[", c(1,3))

chir1 <- read.csv("FilteredLocationCHIR1.0.txt", sep="\t", row.names=1)
chir1 <- cbind(chir1, Pos = (chir1[,"Start"] + chir1[,"Stop"])/2)
chir2 <- read.csv("FilteredLocationCHIR2.0.txt", sep="\t", row.names=1)
chir2 <- cbind(chir2, Pos = (chir2[,"Start"] + chir2[,"Stop"])/2)

snpinfo <- cbind(snpinfo, Chr_C1 = NA)
snpinfo <- cbind(snpinfo, Pos_C1 = NA)

snpinfo <- cbind(snpinfo, Chr_C2 = NA)
snpinfo <- cbind(snpinfo, Pos_C2 = NA)

snpinfo[rownames(chir1), "Chr_C1"] <- chir1[,"chrN"]
snpinfo[rownames(chir1), "Pos_C1"] <- chir1[,"Pos"]

snpinfo[rownames(chir2), "Chr_C2"] <- chir2[,"chrN"]
snpinfo[rownames(chir2), "Pos_C2"] <- chir2[,"Pos"]

if(!file.exists("filtered_snps_numeric.txt")){
  numsnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames = list(rownames(snpdata), colnames(snpdata)))
  for(x in 1:length(snpAlleles)) {
    if(!is.na(snpinfo[x, "reference"]) && snpAlleles[[x]][1] !=  snpinfo[x, "reference"]){  # C/T while reference is T, so flip it around
      snpAlleles[[x]] <- snpAlleles[[x]][2:1]
    }

    g1 <- paste(snpAlleles[[x]][1], snpAlleles[[x]][1],sep="")
    g2a <- paste(snpAlleles[[x]][1], snpAlleles[[x]][2],sep="")
    g2b <- paste(snpAlleles[[x]][2], snpAlleles[[x]][1],sep="")
    g3 <- paste(snpAlleles[[x]][2], snpAlleles[[x]][2],sep="")
    if(!all(snpdata[x,] %in% c(g1,g2a,g2b,g3, NA))) stop("Nope")
    numsnpdata[x, which(snpdata[x, ] == g1)] <- 1
    numsnpdata[x, which(snpdata[x, ] == g2a)] <- 2
    numsnpdata[x, which(snpdata[x, ] == g2b)] <- 2
    numsnpdata[x, which(snpdata[x, ] == g3)] <- 3
  }

  write.table(numsnpdata, "filtered_snps_numeric.txt", sep="\t", quote=FALSE)
}else{
  numsnpdata <- read.csv("filtered_snps_numeric.txt", sep="\t", check.names=FALSE)
}

if(!file.exists("filtered_snps_AB.txt")){
  absnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames = list(rownames(snpdata), colnames(snpdata)))
  for(x in 1:length(snpAlleles)) {
    if(!is.na(snpinfo[x, "reference"]) && snpAlleles[[x]][1] !=  snpinfo[x, "reference"]){  # C/T while reference is T, so flip it around
      snpAlleles[[x]] <- snpAlleles[[x]][2:1]
    }

    g1 <- paste(snpAlleles[[x]][1], snpAlleles[[x]][1],sep="")
    g2a <- paste(snpAlleles[[x]][1], snpAlleles[[x]][2],sep="")
    g2b <- paste(snpAlleles[[x]][2], snpAlleles[[x]][1],sep="")
    g3 <- paste(snpAlleles[[x]][2], snpAlleles[[x]][2],sep="")
    if(!all(snpdata[x,] %in% c(g1,g2a,g2b,g3, NA))) stop("Nope")
    absnpdata[x, which(snpdata[x, ] == g1)]  <- "AA"
    absnpdata[x, which(snpdata[x, ] == g2a)] <- "AB"
    absnpdata[x, which(snpdata[x, ] == g2b)] <- "AB"
    absnpdata[x, which(snpdata[x, ] == g3)]  <- "BB"
  }

  write.table(absnpdata, "filtered_snps_AB.txt", sep="\t", quote=FALSE)
}else{
  absnpdata <- read.csv("filtered_snps_AB.txt", sep="\t", check.names=FALSE)
}

#colnames(numsnpdata) <- paste0(samples[colnames(numsnpdata), "locationShort"], "_", samples[colnames(numsnpdata), "Breed"])

# Add the reference animal to the dataset, since reference is always coded as 1
#numsnpdata <- cbind(numsnpdata, Reference = 1)
clustering <- hclust(dist(t(numsnpdata), method = "manhattan"))

plot(as.phylo(clustering), type="r")  # Or a rooted dendrogram plot: plot(root(as.phylo(clustering),"Reference"), type="r")
plot(clustering, hang = -1)

breeds <- as.character(unique(samples[,"Breed"]))
# Minor Allele Frequencies for different breeds
MAFs <- matrix(NA, nrow(numsnpdata),length(breeds), dimnames=list(rownames(numsnpdata), breeds))
for(breed in breeds){
  individuals <-  rownames(samples)[which(samples[,"Breed"] == breed)]
  MAFs[, breed] <- apply(numsnpdata[, individuals], 1, function(x){
    tabulated <- table(unlist(x))
    ref <- tabulated["1"] * 2 + tabulated["2"]
    alt <- tabulated["3"] * 2 + tabulated["2"]
    if(is.na(ref) || is.na(alt)) return(0)
    if(ref < alt) return(ref / (ref+alt))
    if(ref >= alt) return(alt / (ref+alt))
  })
}

breedSpecific <- names(which(apply(MAFs,1,function(x){(length(which(x == 0)) == 3)})))

## STRUCTURE
if(!file.exists("cleaned_genotypes_structure_NO.txt")){
  numsnpdata <- numsnpdata[,-which(colnames(numsnpdata) == "DN 2")] # Throw away the duplicate individual because it confuses STRUCTURE

  # Write out the data for STRUCTURE
  numsnpdata <- t(numsnpdata)
  structGeno <- NULL #matrix(NA, nrow(numGeno) * 2, ncol(numGeno))
  for(x in 1:nrow(numsnpdata)){
    gg <- rbind(rep(NA, ncol(numsnpdata)), rep(NA, ncol(numsnpdata)))
    a1 <- which(numsnpdata[x,] == 1)
    a2 <- which(numsnpdata[x,] == 2)
    a3 <- which(numsnpdata[x,] == 3)
    gg[1, a1] <- 0; gg[2, a1] <- 0  # Pretty inefficient, but it will do the trick
    gg[1, a2] <- 0; gg[2, a2] <- 1
    gg[1, a3] <- 1; gg[2, a3] <- 1
    gg[is.na(gg)] <- 9
    structGeno <- rbind(structGeno, gg)
  }

  rownames(structGeno) <- gsub(" ","", unlist(lapply(rownames(numsnpdata), rep, 2)))    # Spaces are not allowed in sample names
  colnames(structGeno) <- colnames(numsnpdata)
  write.table(structGeno, file="cleaned_genotypes_structure_NO_DN2.txt", sep = "\t")           # Save the genotypes to disk
}


library(StAMPP)

stammpinput <- t(absnpdata)
stammpinput <- cbind(Sample = rownames(stammpinput), Pop = as.character(samples[rownames(stammpinput),"Breed"]), Ploidy = 2, Format = "BiA", stammpinput)
stammpinput <- as.data.frame(stammpinput)

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values
stammpinput.fst <- stamppFst(stammpinput.freq, 1000, 95, 4) # Population Fst values
stammpinput.fst$Fsts

write.table(stammpinput.fst$Fsts, file = "fsts.txt", sep = "\t")
