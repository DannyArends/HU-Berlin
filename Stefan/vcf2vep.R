# vcf2vep.R - Convert VCF to VEP
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Oct, 2014
# first written Oct, 2014

VCF2VEP <- function(filename = "input.vcf", output = "output.vep", minQUAL = 20, minREADS = 2){
  cat("Reading VCF file\n")
  vcffile <- read.table(filename, colClasses="character")
  colnames(vcffile) <- c("Chr","Start","DBsnp","REF","ALT","QUAL","Filter","INFO","DESC","READS")

  cat("Filtering for low quality reads\n")
  lowQual <- which(as.numeric(vcffile[,"QUAL"]) < minQUAL)
  if(length(lowQual) > 0) vcffile <- vcffile[-lowQual, ]                           # Filter for quality

  cat("Filtering for low read depth\n")
  lowCounts <- which(as.numeric(unlist(lapply(vcffile[,"INFO"], function(x){ 
    splitted <- unlist(strsplit(x,";"))
    substring(splitted[grep("DP", splitted)], 4)
  }))) < minREADS)
  if(length(lowCounts) > 0) vcffile <- vcffile[-lowCounts,]                         # Filter for counts

  cat("Removing complex alternative alleles\n")
  complexALT <- grep(",", vcffile[,"ALT"])                                          # Filter for complex ALT alleles
  if(length(complexALT) > 0) vcffile <- vcffile[-complexALT,]

  cat("Determining SNPs, INSERTIONS and DELETIONS\n")
  snps <- apply(vcffile, 1, function(x){
    if(nchar(x["REF"]) == 1 && nchar(x["ALT"]) == 1) return(TRUE);
    return(FALSE)
  })

  cat("Conversion to VEP format\n")
  vepsnps <- matrix(apply(vcffile[which(snps), c("Chr","Start","REF","ALT")], 1, function(x){
    return(c(x["Chr"],x["Start"],x["Start"],paste0(x["REF"],"/",x["ALT"]),"+"))
  }), length(which(snps)), 5, byrow = TRUE)

  vepother <- matrix(unlist(apply(vcffile[which(!snps), c("Chr","Start","REF","ALT")], 1, function(x){
    startL <- as.numeric(x["Start"])
    if(nchar(x["REF"]) < nchar(x["ALT"])) return(c(x["Chr"], startL+1, startL, paste0("-/", substring(x["ALT"],2)),"+"))                                  # INSERTION (Weird end-position)
    if(nchar(x["REF"]) > nchar(x["ALT"])) return(c(x["Chr"], startL+1, startL + (nchar(x["REF"])-1), paste0(substring(x["REF"],2),"/-"),"+"))             # DELETION
  })), length(which(!snps)), 5, byrow = TRUE)

  cat("Writing the VEP output\n")
  write.table(rbind(vepother, vepsnps), output, sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
}

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
VCF2VEP("Analysis/5073_TTAGGC_L006_.snps.vcf")