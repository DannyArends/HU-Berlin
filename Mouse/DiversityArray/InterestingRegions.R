# InterestingRegions.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of candidate regions from Sebastiaan

setwd("E:/Mouse/DNA/DiversityArray/")

snpOUT  <- read.table("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
snps861S1  <- read.table("Analysis/Diabetes/BFMI861-S1vsALL_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
snpsS1S2  <- read.table("Analysis/Diabetes/BFMI861-S1vsBFMI861-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))


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

regions <- findRegions(snpOUT)
write.table(regions, "Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_Regions.txt", sep="\t", row.names=FALSE)
regions861S1 <- findRegions(snps861S1)
write.table(regions861S1, "Analysis/Diabetes/BFMI861-S1vsALL_Regions.txt", sep="\t", row.names=FALSE)

regionsS1S2 <- findRegions(snpsS1S2)
write.table(regionsS1S2, "Analysis/Diabetes/BFMI861-S1vsBFMI861-S2_Regions.txt", sep="\t", row.names=FALSE)


regionsD1 <- findRegions(deike1)
write.table(regionsD1, "Analysis/Diabetes/Deike1_Regions.txt", sep="\t", row.names=FALSE)

regionsD3 <- findRegions(deike3)
write.table(regionsD3, "Analysis/Diabetes/Deike3_Regions.txt", sep="\t", row.names=FALSE)

# Which markers that we have in the LAB are in the regions we defined
findMarkersInRegions <- function(regions){
  markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
  markersInRegions <- NULL
  for(r in 1:nrow(regions)){
    region <- regions[r,]
    inregion <- which(region["Chr"] == markers[,"Chr"] & region["Start"]-2500000 < markers[,"Location"] & region["End"]+2500000  > markers[,"Location"])
    markersInRegions <- c(markersInRegions, paste(unlist(markers[inregion,"markerID"]),collapse=", "))
  }
  return(cbind(regions,markersInRegions))
}

write.table(findMarkersInRegions(regions), file="Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_MarkersInRegions.txt", sep="\t", row.names=FALSE)
write.table(findMarkersInRegions(regions861S1), file="Analysis/Diabetes/BFMI861-S1vsALL_MarkersInRegions.txt", sep="\t", row.names=FALSE)

query.regions <- paste(regions[,"Chr"],regions[,"Start"],regions[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]


write.table(res.biomart, "Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_GenesInRegions.txt", sep="\t", row.names=FALSE)