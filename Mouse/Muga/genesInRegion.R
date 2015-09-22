
genesInRegion <- function(chr, start_position, end_position, species = "mmusculus", verbose = TRUE){
  require("biomaRt")
  if(is.na(end_position)) end_position <- start_position + 1
  mart = useMart("ensembl", dataset = paste0(species, "_gene_ensembl"))
  chr.region = c(paste0(chr, ":",start_position, ":", end_position))
  filterlist = list(chr.region, "protein_coding")
  if(verbose) cat("Retrieving results for ", chr.region, "\n")
  res <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"), filters = c("chromosomal_region","biotype"), values = filterlist, mart = mart)
  if(verbose) cat("Retrieved ",nrow(res), "results for ", chr.region, "\n")
  return(res)
}

setwd("E:/Mouse/DNA/MegaMuga/")

patRegions <- read.table("Analysis/patRegionsBias.txt",sep="\t",header=TRUE, check.names=FALSE)
results <- NULL
for(x in 1:nrow(patRegions)){
    results <- rbind(results, genesInRegion(patRegions[x,"Chr"],patRegions[x,"Start"],patRegions[x,"Stop"]))
}
cat(unique(results[,"ensembl_gene_id"]),sep="\n", file="Analysis/genesPatBias.txt")
#Retrieved  12 results for  2:3164247:4280400 
#Retrieved  31 results for  3:6274425:14841429 
#Retrieved  68 results for  4:3569913:20188902 
#Retrieved  27 results for  4:31935558:35049506 
#Retrieved  3 results for  4:75974884:81245565 
#Retrieved  182 results for  4:112951834:125311291 
#Retrieved  38 results for  5:3160082:9746133 
#Retrieved  1 results for  5:127617133:127617134 
#Retrieved  53 results for  6:76542612:84179155 
#Retrieved  8 results for  6:94268595:96800636 
#Retrieved  57 results for  9:75916734:86816288 
#Retrieved  12 results for  10:9084536:12411896 
#Retrieved  16 results for  12:5253913:10714469 
#Retrieved  32 results for  14:7705302:14880907 
#Retrieved  1 results for  14:92378135:93875670 
#Retrieved  1 results for  16:5617528:7969824 
#Retrieved  37 results for  16:10390703:14151479 
#Retrieved  1 results for  16:97452473:97452474 
#Retrieved  3 results for  17:5991544:6097036 
#Retrieved  112 results for  17:33612453:35503627 
#Retrieved  35 results for  17:47490686:50134302 
#Retrieved  20 results for  18:10553145:13902153 
#Retrieved  23 results for  18:38281545:42987483 
#Retrieved  17 results for  19:10097832:10662901 
#Retrieved  4 results for  19:24458355:24741711 

matRegions <- read.table("Analysis/matRegionsBias.txt",sep="\t",header=TRUE, check.names=FALSE)
results <- NULL
for(x in 1:nrow(matRegions)){
    results <- rbind(results, genesInRegion(matRegions[x,"Chr"],matRegions[x,"Start"],matRegions[x,"Stop"]))
}
cat(unique(results[,"ensembl_gene_id"]),sep="\n", file="Analysis/genesMatBias.txt")
#Retrieved  12 results for  1:3668628:6357478 
#Retrieved  12 results for  2:3164247:4280400 
#Retrieved  4 results for  2:55092420:57658939 
#Retrieved  31 results for  3:6274425:14841429 
#Retrieved  26 results for  4:3569913:11002647 
#Retrieved  54 results for  5:3160082:13395581 
#Retrieved  15 results for  6:3371606:5248452 
#Retrieved  4 results for  6:21943927:22463199 
#Retrieved  15 results for  6:96800636:100429064 
#Retrieved  21 results for  9:74576994:77181656 
#Retrieved  2 results for  9:100279699:100523422 
#Retrieved  16 results for  10:12091154:15704536 
#Retrieved  3 results for  11:29036920:29318982 
#Retrieved  14 results for  12:5253913:9866543 
#Retrieved  1 results for  12:100937187:100937188 
#Retrieved  76 results for  13:17868992:22380465 
#Retrieved  14 results for  14:7705302:9716964 
#Retrieved  122 results for  16:5617528:18450764 
#Retrieved  43 results for  17:3737834:11318508 
#Retrieved  1 results for  17:33661748:33661749 
#Retrieved  44 results for  18:4516519:13902153 
#Retrieved  17 results for  19:22845930:24741711 
