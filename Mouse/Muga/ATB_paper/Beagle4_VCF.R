#
# Analyse megaMuga data from multiple generations using beagle and R
# Phase genotype data using BEAGLE 4.0
#

source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")
source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/vcfTools.R")

# Load the genotype call data
setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
locusxdnaheader <- unlist(strsplit(readLines("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", n=16)[16],","))
locusxdnaR <- read.csv("Humboldt_Berlin_MEGMUGV01_20140817/Humboldt_Berlin_MEGMUGV01_20140817_LocusXDNA.csv", header=FALSE, skip = 22)
colnames(locusxdnaR) <- c("Label","plateWell","Date","oligoPoolId","bundleId", "status", "Type", "Nas", locusxdnaheader[4:length(locusxdnaheader)])
locusxdna <- t(locusxdnaR[which(locusxdnaR[,"Type"] == "calls"),])
locusxdnascores <- t(locusxdnaR[which(locusxdnaR[,"Type"] == "Score_Call"),])
colnames(locusxdna) <- locusxdna[1,]
colnames(locusxdnascores) <- locusxdnascores[1,]
locusxdna <- locusxdna[-c(1:8), ]                                                     # Remove some unneeded headers
locusxdnascores <- locusxdnascores[-c(1:8), ]                                         # Remove some unneeded headers
locusxdna[which(locusxdna == "U")] <- NA                                              # Use the correct NA values
#locusxdna[which(locusxdnascores < 0.7)] <- NA                                         # Set NA unreliable calls
# Remove markers and individuals without any data and repeated measurements for the same marker
locusxdna <- locusxdna[-which(apply(locusxdna,1,function(x){return(sum(is.na(x)))}) == ncol(locusxdna)),]
locusxdna <- locusxdna[,-which(apply(locusxdna,2,function(x){return(sum(is.na(x)))}) == nrow(locusxdna))]

cat("Loaded:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 479 individuals, 59781 markers
if(!file.exists("TRDredo/AllGenotypes.txt")) write.table(locusxdna, "TRDredo/AllGenotypes.txt", sep="\t", quote=FALSE)

samples <- colnames(locusxdna)
markers <- rownames(locusxdna)

#Load annotation from JAX
if(!file.exists("TRDredo/markerAnnotation.txt")){
  load("MM_snps.Rdata")
  MM_snps[,"SNP_ID"] <- gsub(".", "-", MM_snps[,"SNP_ID"], fixed = TRUE)                # Fix the . to - error in the MM_SNP
  rownames(MM_snps) <- MM_snps[,"SNP_ID"]                                               # Use them as rownames

  marker.annot <- matrix(NA, length(markers), 4, dimnames=list(markers,c("Chr", "Pos", "snp_A", "snp_B")))
  marker.annot[markers, "Chr"] <- MM_snps[markers,"Chr"]
  marker.annot[markers, "Pos"] <- MM_snps[markers,"Mb_NCBI38"] * 1000000

  for(x in markers){
    instr  <- regexpr("\\[.+\\]", MM_snps[x,"Sequence"], perl=TRUE)
    marker.annot[x,c("snp_A","snp_B")] <- unlist(strsplit(substr(MM_snps[x,"Sequence"], (instr[1]+1), (instr[1]+attr(instr,"match.length")-2)),"/"))
  }
  write.table(marker.annot, "TRDredo/markerAnnotation.txt", sep="\t", quote=FALSE)
}else{
  marker.annot <- read.table("TRDredo/markerAnnotation.txt", sep="\t")
}

# Remove individuals with wrong parentals
locusxdna <- locusxdna[,-which(colnames(locusxdna) %in% c("6661965", "6662155", "6662156","6661257", "6661258", "6661259", "6661260", "6661261", "33310233"))]

excludedM <- c()

# Remove the duplicated markers
duplicates <- rownames(marker.annot)[which(duplicated(apply(marker.annot[,c("Chr","Pos")],1,paste0,collapse="")))]
excludedM <- cbind(marker.annot[which(rownames(marker.annot) %in% duplicates),], reason = "duplicate")
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% duplicates),]
locusxdna <- locusxdna[-which(rownames(locusxdna) %in% duplicates),]

cat("After duplicate removal:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 470 individuals, 70597 markers

# Remove markers that do not have a known location
missingloc <- rownames(marker.annot)[which(is.na(marker.annot[,"Pos"]))]
excludedM <- rbind(excludedM, cbind(marker.annot[which(rownames(marker.annot) %in% missingloc),], reason = "noposition"))
if(length(missingloc) > 0){
  marker.annot <- marker.annot[-which(rownames(marker.annot) %in% missingloc),]
  locusxdna <- locusxdna[-which(rownames(locusxdna) %in% missingloc),]
}

cat("After duplicate and unplaced marker removal:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 470 individuals, 70594 markers

missingPerMarker <- apply(locusxdna, 1, function(x){ sum(is.na(x)) / length(x) * 100 })
tooMuchMissing <- names(which(missingPerMarker > 10))
excludedM <- rbind(excludedM, cbind(marker.annot[which(rownames(marker.annot) %in% tooMuchMissing),], reason = "missing"))
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% tooMuchMissing),]
locusxdna <- locusxdna[-which(rownames(locusxdna) %in% tooMuchMissing),]
cat("After duplicate and unplaced marker and missing removal:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 470 individuals, 70390 markers

# Remove markers where B6N and BFMI aren't homozygous
wrongOrigin <- rownames(locusxdna)[which(is.na(locusxdna[,"B6N"]) | is.na(locusxdna[,"BFMI860-12 (V2)"]) | locusxdna[,"B6N"] =="H" | locusxdna[,"BFMI860-12 (V2)"] =="H")]
excludedM <- rbind(excludedM, cbind(marker.annot[which(rownames(marker.annot) %in% wrongOrigin),], reason = "founderequal"))
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% wrongOrigin),]
locusxdna <- locusxdna[-which(rownames(locusxdna) %in% wrongOrigin),]

cat("After duplicate and unplaced marker and missing and founder origin removal:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 470 individuals, 65950 markers

# Remove markers where B6N == BFMI 
founderEqual <- rownames(locusxdna)[which(locusxdna[,"B6N"] == locusxdna[,"BFMI860-12 (V2)"])]
excludedM <- rbind(excludedM, cbind(marker.annot[which(rownames(marker.annot) %in% founderEqual),], reason = "founderequal"))
marker.annot <- marker.annot[-which(rownames(marker.annot) %in% founderEqual),]
locusxdna <- locusxdna[-which(rownames(locusxdna) %in% founderEqual),]
cat("After duplicate and unplaced marker and missing and founder origin and equality removal:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 470 individuals, 18329 markers

## Order the markers since the VCF needs to have continuous markers per chromosome
marker.annot <- marker.annot[order(marker.annot[,"Chr"], as.numeric(marker.annot[,"Pos"])),]
markers <- rownames(marker.annot)
locusxdna <- locusxdna[markers, ]
cat("Final data set for phasing:", ncol(locusxdna), "individuals,", nrow(locusxdna), "markers\n")  ###Loaded: 470 individuals, 18329 markers

write.table(marker.annot, "TRDredo/markers.annot", sep="\t", quote=FALSE)
write.table(excludedM, "TRDredo/marker.excluded", sep="\t", quote=FALSE)
write.table(locusxdna, "TRDredo/genotypes.txt", sep="\t", quote=FALSE)


samples <- colnames(locusxdna)
markers <- rownames(locusxdna)

# Write genotypes as VCF for Beagle4
if(!file.exists("TRDredo/genotypes.vcf")){
  cat("##fileformat=VCFv4.2\n", file="TRDredo/genotypes.vcf")
  cat("##FORMAT=<ID=GT,Number=1,Type=Integer,Description=\"Genotype\">\n", file="TRDredo/genotypes.vcf", append = TRUE)
  cat(paste0("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t", paste(samples, collapse = "\t"), "\n"), file="TRDredo/genotypes.vcf", append = TRUE)
  for(x in markers){
    gts <- paste(toVCFcode(locusxdna[x, samples]), collapse = "\t")
    cat(marker.annot[x,"Chr"], "\t", marker.annot[x,"Pos"], "\t", x, "\t", marker.annot[x,"snp_A"], "\t", marker.annot[x,"snp_B"], "\t.\tPASS\t.\tGT\t", gts, "\n", sep="", file="TRDredo/genotypes.vcf", append = TRUE)
  }
}

# Load in the phenotypes and create the pedigree file
phenotypes <- read.table("Phenotypes/allPhenotypes.txt",sep="\t",header=TRUE,na.strings=c(0,"-","NA"))
rownames(phenotypes) <- phenotypes[,"ID"]
phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix

if(!file.exists("TRDredo/population.ped")){
  cat("", file="TRDredo/population.ped") # Empty ped file holding the pedigree
  for(x in samples){
    dataline <- c("1", x, phenotypes[x,"Vater"], phenotypes[x,"Mutter"])
    dataline[is.na(dataline)] <- "0"
    cat(paste(dataline, collapse=" "), "\n", sep = "", file="TRDredo/population.ped", append=TRUE) # Empty ped file holding the pedigree
  }
}

if(!file.exists("TRDredo/phased.vcf.gz")){ # Do beagle phasing
  system(paste0("java -Xmx1000m -jar beagle.v.4.jar gt=TRDredo/genotypes.vcf ped=TRDredo/population.ped out=TRDredo/phased"))
}
