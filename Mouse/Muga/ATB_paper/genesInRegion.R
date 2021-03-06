genesInRegion <- function(chr, start_position, end_position, species = "mmusculus", verbose = TRUE){
  require("biomaRt")
  if(is.na(end_position)) end_position <- start_position + 1
  mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="http://nov2020.archive.ensembl.org")
  chr.region = c(paste0(chr, ":",start_position, ":", end_position))
  filterlist = list(chr.region, "protein_coding")
  if(verbose) cat("Retrieving results for ", chr.region, "\n")
  res <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "mgi_symbol", "mgi_description"), filters = c("chromosomal_region","biotype"), values = filterlist, mart = mart)
  if(verbose) cat("Retrieved ",nrow(res), "results for ", chr.region, "\n")
  return(res)
}

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")

# Load the Allele Transmission Biased regions
patRegions <- read.table("TRDredo/pat0.01annot.txt", sep="\t", header=FALSE)
matRegions <- read.table("TRDredo/mat0.01annot.txt", sep="\t", header=FALSE)
colnames(patRegions) <- c("Code", "Flanking.Marker", "Chr", "Start", "TopPos", "Stop",  "Flanking.Marker.1", "Size", "Prefered Allele", "nSNPs", "T1", "T2", "Top", "TopT1", "TopT2", "perc")
colnames(matRegions) <- c("Code", "Flanking.Marker", "Chr", "Start", "TopPos", "Stop",  "Flanking.Marker.1", "Size", "Prefered Allele", "nSNPs", "T1", "T2", "Top", "TopT1", "TopT2", "perc")

BFMIregionsPAT <- patRegions[which(patRegions[,"Prefered Allele"] == "BFMI"),]
B6NregionsPAT <- patRegions[which(patRegions[,"Prefered Allele"] == "B6N"),]

BFMIregionsMAT <- matRegions[which(matRegions[,"Prefered Allele"] == "BFMI"),]
B6NregionsMAT <- matRegions[which(matRegions[,"Prefered Allele"] == "B6N"),]

doSearch <- function(regions, filestem){
  results <- NULL
  for(x in 1:nrow(regions)){
    results <- rbind(results, genesInRegion(regions[x,"Chr"], regions[x,"Start"], regions[x,"Stop"]))
  }
  cat(unique(results[,"ensembl_gene_id"]),sep="\n", file=paste0("TRDredo/",filestem,"genes.txt"))
  write.table(results, file=paste0("TRDredo/",filestem, "genesInfo.txt"),sep="\t", quote = FALSE, row.names = FALSE)
}

doSearch(BFMIregionsPAT, "PAT_BFMI")
doSearch(B6NregionsPAT, "PAT_B6N")
doSearch(BFMIregionsMAT, "MAT_BFMI")
doSearch(B6NregionsMAT, "MAT_B6N")
