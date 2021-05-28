
# retrieve entrez gene ID
require("biomaRt")
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl", host="http://nov2020.archive.ensembl.org")
ensembltoentrez <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = c("ensembl_gene_id"), values = unique(as.character(allgenes[,1])), mart = mart)
write.table(ensembltoentrez, "TRDredo/genesToEntrez.txt",sep="\t")

# Look at P:P interactions in BIOGRID
biogriddata <- read.csv("stringdb/BIOGRID-ALL-3.4.132.tab2.txt",sep="\t")
biogriddata <- biogriddata[which(biogriddata[,"Organism.Interactor.A"] == 10090 & biogriddata[,"Organism.Interactor.B"] == 10090),]
biogriddata <- biogriddata[which(as.character(biogriddata[,"Entrez.Gene.Interactor.A"]) %in% as.character(ensembltoentrez[,2]) & biogriddata[,"Entrez.Gene.Interactor.B"] %in% ensembltoentrez[,2]),]

# Create results
interacting <- biogriddata[which(as.character(biogriddata[,c(2,3,8,9)][,3]) != as.character(biogriddata[,c(2,3,8,9)][,4])), c(2,3,8,9)]
interacting <- cbind(interacting, RegionGeneA = NA, RegionGeneB = NA, ChiSq = NA, LOD = NA)

# Which region the gene is in
geneLocation <- function(allgenes, ensgid, verbose = TRUE){
  ii <- allgenes[which(allgenes[,1] == ensgid),]
  cat(ensgid, ": Found", nrow(ii), "matching genes\n")
  gChr <- ii[1,"chromosome_name"];gLoc <- ii[1,"start_position"]; gLoc1 <- ii[1,"end_position"]  # Chr:Location
  cat(ensgid, ": ", gChr, "-", gLoc, "\n")
  return(paste0(gChr,":", round(gLoc / 1000000,d = 2),"-", round(gLoc1/ 1000000, d = 2)))
}


# Which region the gene is in
whichRegion <- function(allgenes, ensgid, verbose = TRUE){
  ii <- allgenes[which(allgenes[,1] == ensgid),]
  cat(ensgid, ": Found", nrow(ii), "matching genes\n")
  gChr <- ii[1,"chromosome_name"];gLoc <- ii[1,"start_position"]  # Chr:Location
  cat(ensgid, ": ", gChr, "-", gLoc, "\n")
  rownames(allRegions[which(allRegions[,"Chr"] == gChr & as.numeric(allRegions[,"Start"]) < gLoc & as.numeric(allRegions[,"Stop"]) > gLoc ),])
}

## Get regions and Chi2 scores
for(x in 1:nrow(interacting)){
  ensgidA <- ensembltoentrez[which(ensembltoentrez[,2] == interacting[x,1]),1]
  ensgidB <- ensembltoentrez[which(ensembltoentrez[,2] == interacting[x,2]),1]
  interacting[x, "RegionGeneA"] <- whichRegion(allgenes, ensgidA)[1]
  interacting[x, "RegionGeneB"] <- whichRegion(allgenes, ensgidB)[1]
  if(!is.na(interacting[x, "RegionGeneA"]) && !is.na(interacting[x, "RegionGeneB"])){
    interacting[x, "ChiSq"] <- results[ interacting[x, "RegionGeneA"] , interacting[x, "RegionGeneB"] ]
    interacting[x, "LOD"] <- LODscores[ interacting[x, "RegionGeneA"] , interacting[x, "RegionGeneB"] ]
  }
}
colnames(interacting)[1:4] <- c("Entrez.A","Entrez.B","Symbol.A", "Symbol.B")

