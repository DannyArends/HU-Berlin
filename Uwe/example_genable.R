# Use Siham's Goat data and load it into GenABEL for QTL scanning, an example for Uwe
#
# copyright (c) 2016-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jul, 2016
# first written Jul, 2016

# Uncomment the following line if you do not have genable installed in R
#install.packages("GenABEL")
library(GenABEL)

# Set the working directory to where the data is stored, replace \ by / in windows
setwd("E:/UWE")

# Load the raw SNP data
snpdata <- read.csv("Ziegen_HU-Berlin_Matrix.txt", sep="\t", skip = 9, header=TRUE, check.names=FALSE, row.names=1, na.strings=c("NA", "", "--"))
# SNP information (Allele, locations etc)
snpinfo <- read.table("snpinfo.txt", sep="\t", na.strings=c("", "NA", "N.D."))

snpdata <- snpdata[rownames(snpinfo), ]                         # Take only the SNPs for which we have information (2000 here)
snpdata <- snpdata[,-which(colnames(snpdata) == "DN 2")]        # Throw away the duplicate individual because it confuses STRUCTURE

# Load the phenotype data for the samples
samples    <- read.table("sampleinfo.txt", sep="\t")
# Load the fixed location data
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`
samples    <- cbind(samples, locations[rownames(samples),])
samples    <- cbind(samples, locationShort = as.character(unlist(lapply(strsplit(as.character(samples[,"Location"]), "_"), "[",1))))
rownames(samples) <- gsub(" ", "", rownames(samples))                         # Sample names cannot contain spaces

gendata <- cbind(snpinfo[rownames(snpdata),c("Chr", "Position")], snpdata)    # Add location information to the genotype data
colnames(gendata) <- gsub(" ", "", colnames(gendata))                         # Sample names cannot contain spaces

# Create the phenotype and covariate file we want to use: sex and age are required ? I add averagemilk as phenotype for QTL scanning
phenocovs <- cbind(id = rownames(samples), sex = rep(0, nrow(samples)), age = samples[,"Age"], averagemilk = samples[,"Averagemilk"])

# Write the genotypes
write.table(cbind(id = rownames(gendata), gendata), file="genable.input", sep="\t", quote=FALSE, na="00", row.names=FALSE)
# Write the phenotypes
write.table(phenocovs, file="genable.pheno", sep="\t", quote=FALSE)

# Convert them to genable binary format using the provided function for csv
convert.snp.illumina("genable.input", "genable.encoded", strand = "+")

# Load the converted genotype data, and the phenotype data
gwadata <- load.gwaa.data(phenofile = "genable.pheno", genofile = "genable.encoded", force = TRUE, makemap = FALSE, sort = TRUE, id = "id")

# Scan for a QTL, adjusted for sex and age (CRSNP is the currentSNP in the model)
res <- scan.glm("averagemilk ~ sex + age + CRSNP", data=gwadata)
