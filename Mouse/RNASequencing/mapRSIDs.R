# Map RS id's to the ASE data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015
#

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/Paper")
mdata <- read.table("S3_ASE.txt", sep="\t", header=TRUE)
positions <- paste0(mdata[,"CHROM"],":",mdata[,"POS"],":",mdata[,"POS"])          # Use the positions to request SNPs from biomart
snps <- read.table("results.txt", sep="\t", header=TRUE)                          # Obtained from biomart (web)

mdata <- cbind(mdata, snpinfo = NA)
for(x in 1:nrow(mdata)){
  subsnps <- snps[which(snps[,3] == mdata[x,"CHROM"]),]
  rsid <- as.character(subsnps[subsnps[,2] ==  mdata[x,"POS"],1])
  if(length(rsid) > 0){
  mdata[x, "snpinfo"] = paste(rsid, collapse="|")
  }
}

write.table(mdata, "S3_ASE_DBsnp.txt",sep="\t", row.names=FALSE)
mdata <- read.table("S3_ASE_DBsnp.txt", sep="\t", header=TRUE, colClasses="character")

vcfdata <- read.table("vep860v2.pes.s.rmd.rg.realigned.recal.haplotclld.vcf", sep = "\t", header=TRUE)

for(x in which(is.na(mdata[,"snpinfo"]))){       # Look if the SNPs are found in our DNA seq data
  chrInVCF <- paste0("chr",mdata[x,"CHROM"])
  subvcf <- vcfdata[which(vcfdata[,"CHROM"] == chrInVCF),]
  hasSNP <- which(as.numeric(subvcf[,"POS"]) == as.numeric(mdata[x,"POS"]))
  if(length(hasSNP) > 0){
    mdata[x,"snpinfo"] <- "Seq"
    cat("Found\n")
  }
}

write.table(mdata, "S3_ASE_DBsnp_VCF.txt",sep="\t", row.names=FALSE)