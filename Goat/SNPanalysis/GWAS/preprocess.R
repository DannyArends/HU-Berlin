#
# Preprocess the goat data for GWAS
#
library(heterozygous)
library(genetics)

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")

#Data is filtered for: GenTrain > 0.6, MAF > 0.05, unknown positions, 5% call rate 
snpdata    <- read.table("filtered_snps_NO_DN2.txt", sep="\t", check.names=FALSE, colClasses="character")
snpinfo    <- read.table("snpinfo.txt", sep="\t", na.strings=c("", "NA", "N.D."))
samples    <- read.table("sampleinfo.txt", sep="\t")
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`

dim(snpdata)

# Filter MAF more aggresively <= 0.1
allelefreq <- apply(snpdata, 1 , function(x){
  tbl <- table(unlist(lapply(x, strsplit, "")))
  min(tbl / sum(tbl))
})

snpinfo <- cbind(snpinfo, MAF=allelefreq)

tooLowAllele <- which(allelefreq <= 0.1)
snpdata <- snpdata[-tooLowAllele, ]

# Filter HWE more aggresively <= 0.1
hwe <- HWE(snpdata)
hweadj <- p.adjust(hwe, "BH")
outofhwe <- which(hweadj < 0.1)

snpdata <- snpdata[-outofhwe, ]

# Scale down the SNPinfo based on the remaining SNP data
snpinfo <- snpinfo[rownames(snpdata), ]

# Filter GenTrain more aggresively > 0.8
lowGenTrain <- which(snpinfo[,"GenTrain.Score"] <= 0.8)
snpinfo <- snpinfo[-lowGenTrain, ]
snpdata <- snpdata[rownames(snpinfo),]

dim(snpdata)
dim(snpinfo)

# Prune MAP for markers in LD R2 > 0.5, keep the ones with the highest MAF
map <- snpinfo[,c("Chr", "Position", "MAF")]
map.pruned <- map

genotypes <- snpdata
genotypes.pruned <- genotypes
markers <- rownames(genotypes.pruned)

x <- 1
while(x < length(markers)) {
  mName <- markers[x]
  mChr <- as.character(map[mName, "Chr"])
  mPos <- as.numeric(map[mName, "Position"])

  nearby <- which(as.character(map.pruned[,"Chr"]) == mChr & as.numeric(map.pruned[, "Position"]) > (mPos - 1000000) &  as.numeric(map.pruned[, "Position"]) < (mPos + 1000000))
  nearby <- rownames(map.pruned)[nearby]

  mGeno <- genotype(as.character(genotypes.pruned[mName,]), sep = "")
  locOfMarker <- which(nearby == mName)
  if(locOfMarker > 1) nearby <- nearby[-(1:locOfMarker-1)]

  LDs <- rep(NA, length(nearby))
  names(LDs) <- nearby
  for(y in nearby){
     LDs[y] <- round(LD(mGeno, genotype(as.character(genotypes.pruned[y,]), sep = ""))$"R^2",2)
  }
  inLD <- names(which(LDs > 0.5))
  if(length(inLD) > 1){
    MAFs <- map.pruned[inLD, ]
    bestSNP <- rownames(MAFs)[which.max(MAFs[, "MAF"])]
    inLD <- inLD[-which(inLD == bestSNP)]
    toPrune <- which(rownames(genotypes.pruned) %in% inLD)
    if(length(toPrune) > 0){
      genotypes.pruned <- genotypes.pruned[-toPrune, ]
      map.pruned <- map.pruned[-toPrune, ]
      markers <- rownames(genotypes.pruned)
    }
  }
  cat("Done", x, "/", nrow(genotypes), "/", nrow(genotypes.pruned), "==", nrow(map.pruned), "\n")
  x <- (x + 1)
}

# Fix up the samples, and add the phenotypes from the location object
samples <- samples[colnames(genotypes.pruned),]
samples <- cbind(samples, locations[rownames(samples),])

# We demand at least 10 individuals in each genotype class
NminAllele <- unlist(lapply(apply(genotypes.pruned,1,table),min))
genotypes.pruned <- genotypes.pruned[-which(NminAllele < 10),]
map.pruned <- map.pruned[rownames(genotypes.pruned),]

# Write out the files for GWAS
setwd("D:/Edrive/Goat/DNA/SihamAnalysis/GWAS")
write.table(genotypes.pruned, "genotypes.txt", sep="\t", quote=FALSE)
write.table(map.pruned, "geneticmap.txt", sep="\t", quote=FALSE)
write.table(samples, "samples.txt", sep="\t", quote=FALSE)

