# Analysis of Chicken 600K SNP chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

setwd("E:/Chicken/DNA/600KSNPChip/")

arrayAnnotation <- read.table("Annotation/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", header=TRUE, colClasses="character", na.string="---")

library(biomaRt)                                                                                    # Biomart
snp.db <- useMart("snp", dataset="ggallus_snp")                                                     # For chicken SNPs
snps <- as.character(arrayAnnotation[,"dbSNP.RS.ID"])                                               # The list of RS_IDs we want to retrieve
snps <- snps[-which(is.na(snps))]                                                                   # Remove the empty ones

if(!file.exists("Analysis/DBSNPannotation.txt")){
  biomartResults <- NULL
  for(x in seq(1, length(snps), 1000)){                                                             # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(snps))                                                            # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    res.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start"),                          # Use biomart to retrieve locations and reference alleles
                       filters="snp_filter", values=snps[x:xend], mart=snp.db)
    biomartResults <- rbind(biomartResults, res.biomart)
  }
  write.table(biomartResults, file="Analysis/DBSNPannotation.txt", sep="\t", row.names=FALSE)
}else{                                                                                              # If the annotation file is there use it
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("Analysis/DBSNPannotation.txt", sep="\t", header=TRUE)
}
