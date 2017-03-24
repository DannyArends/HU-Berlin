# Deike.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2017
# first written Jan, 2017
#
# Analysis of Diversity array data for Deike

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")

#Deike 1
deike1  <- read.table("Analysis/Diabetes/BFMI861-S2vs861S1n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))

#Deike 2
deike <- read.table(file="Analysis/Diabetes/BFMI861-S2nB6Nvs861S1n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE)


#Deike 3
deike3  <- read.table("Analysis/Diabetes/BFMI861-S1vs861S2n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))



findRegions <- function(snps, snpInRegion = 35, maxDistance = 4000000){
  chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

  snpInRegion <- 35                                                         # A region should contain at least 35 SNPs
  maxDistance <- 4000000                                                    # If the neighbouring SNP is less that 4mb away, we continue with our region

  regions <- NULL
  for(chr in chromosomes){
    onchr  <- snps[snps[,"Chr"] == chr,]
    if(nrow(onchr) > snpInRegion){
      sorted <- sort(as.numeric(onchr[,"Location"]), index.return=TRUE)
      onchr  <- onchr[sorted$ix,]
      sloc   <- 0
      snpcnt <- 0
      for(snp in 1:(nrow(onchr)-1)){
        cloc <- onchr[snp, "Location"]
        nloc <- onchr[snp+1, "Location"]
        if(as.numeric(nloc) - as.numeric(cloc) < maxDistance){
          if(snpcnt == 0) sloc <- cloc
          snpcnt <- snpcnt + 1
        }else{
          if(snpcnt > snpInRegion){
            regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
            cat("New region:", chr, sloc, cloc, snpcnt,"\n")
          }
          snpcnt <- 0
        }
      }
      if(snpcnt > snpInRegion){
        regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
        cat("New region:",chr, sloc, cloc, snpcnt,"\n")
      }
    }
  }

  colnames(regions) <- c("Chr","Start","End","SNPs")
  return(regions)
}

regionsD1 <- findRegions(deike1)
write.table(regionsD1, "Analysis/Diabetes/S2vsAll_Regions.txt", sep="\t", row.names=FALSE)

regionsD3 <- findRegions(deike3)
write.table(regionsD3, "Analysis/Diabetes/S1vsAll_Regions.txt", sep="\t", row.names=FALSE)


regions <- findRegions(deike)
query.regions <- paste(regions[,"Chr"],regions[,"Start"],regions[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]

write.table(res.biomart, "Analysis/Diabetes/BFMI861-S2nB6Nvs861S1n860-12n860-S2_GenesInRegions.txt", sep="\t", row.names=FALSE)
write(res.biomart[,"ensembl_gene_id"], "Analysis/Diabetes/BFMI861-S2nB6Nvs861S1n860-12n860-S2_MGIsymbols.txt")




regions <- findRegions(deike1)

query.regions <- paste(regions[,"Chr"],regions[,"Start"],regions[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]

write.table(res.biomart, "Analysis/Diabetes/BFMI861-S2vsAll_GenesInRegions.txt", sep="\t", row.names=FALSE)


regions <- findRegions(deike3)

query.regions <- paste(regions[,"Chr"],regions[,"Start"],regions[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]

write.table(res.biomart, "Analysis/Diabetes/BFMI861-S1vsAll_GenesInRegions.txt", sep="\t", row.names=FALSE)


