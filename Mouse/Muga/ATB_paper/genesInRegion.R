genesInRegion <- function(chr, start_position, end_position, species = "mmusculus", verbose = TRUE){
  require("biomaRt")
  if(is.na(end_position)) end_position <- start_position + 1
  mart = useMart("ENSEMBL_MART_ENSEMBL")
  dataset = useDataset("mmusculus_gene_ensembl", mart)
  chr.region = c(paste0(chr, ":",start_position, ":", end_position))
  filterlist = list(chr.region, "protein_coding")
  if(verbose) cat("Retrieving results for ", chr.region, "\n")
  res <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "mgi_symbol", "mgi_description"), filters = c("chromosomal_region","biotype"), values = filterlist, mart = dataset)
  if(verbose) cat("Retrieved ",nrow(res), "results for ", chr.region, "\n")
  return(res)
}

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")

patRegions <- read.table("Analysis/ATB_PAT.txt",sep="\t",header=TRUE, check.names=FALSE)
matRegions <- read.table("Analysis/ATB_MAT.txt",sep="\t",header=TRUE, check.names=FALSE)

allRegions <- rbind(patRegions, matRegions)

BFMIregionsPAT <- patRegions[which(patRegions[,"Prefered Allele"] == "BFMI"),]
B6NregionsPAT <- patRegions[which(patRegions[,"Prefered Allele"] == "B6N"),]

BFMIregionsMAT <- matRegions[which(matRegions[,"Prefered Allele"] == "BFMI"),]
B6NregionsMAT <- matRegions[which(matRegions[,"Prefered Allele"] == "B6N"),]

doSearch <- function(regions, filestem){
  results <- NULL
  for(x in 1:nrow(regions)){
    results <- rbind(results, genesInRegion(regions[x,"Chr"], regions[x,"Start"], regions[x,"Stop"]))
  }
  cat(unique(results[,"ensembl_gene_id"]),sep="\n", file=paste0("Analysis/",filestem,"genes.txt"))
  write.table(results, file=paste0("Analysis/",filestem, "genesInfo.txt"),sep="\t", quote = FALSE, row.names = FALSE)
}

doSearch(BFMIregionsPAT, "PAT_BFMI")
doSearch(B6NregionsPAT, "PAT_B6N")
doSearch(BFMIregionsMAT, "MAT_BFMI")
doSearch(B6NregionsMAT, "MAT_B6N")
