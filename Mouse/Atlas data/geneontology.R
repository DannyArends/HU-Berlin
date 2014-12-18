# Analysis of the micro array data from Atlas 2014
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

library(biomaRt)                                                                                        # BiomaRt package
library(topGO)                                                                                          # topGO package

setwd("E:/Mouse/RNA/ArrayDesign/Atlas data")
arrays <- read.table("Annotation/arrays.txt", header=TRUE, sep="\t", colClasses="character")
alldata <- read.table("Analysis/geneexpression.txt", sep="\t", header=TRUE)

DEthreshold <- 0.1 / nrow(alldata)

DEstrain    <- which(as.numeric(alldata[,"Strain_P"]) < DEthreshold)
DEtissue    <- which(as.numeric(alldata[,"Tissue_P"]) < DEthreshold)

cat("Found", length(DEtissue), "probes differentially expressed between tissues\n")
cat("Found", length(DEstrain), "probes differentially expressed between strains\n")

allgenes <- unique(alldata[,"ensembl_gene_id"])

if(!file.exists("Annotation/GOannotation.txt")){
  bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                              # BiomaRt for mouse genes
  biomartResults <- NULL
  for(x in seq(1, length(allgenes), 1000)){                                                             # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(allgenes))                                                            # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    res.biomart <- getBM(c("ensembl_gene_id","go_id"),                                                  # Use biomaRt to retrieve GO terms
                          filters="ensembl_gene_id", values=allgenes[x:xend], mart=bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
  }
  write.table(biomartResults, file="Annotation/GOannotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading biomaRt gene ontology annotation from disk\n")
  biomartResults <- read.table("Annotation/GOannotation.txt", sep="\t", header=TRUE)
}

if(!file.exists("Annotation/geneid2go.map")){                                                         # Create the ENSG to GO map only if it doesn't exists
  cat("", file="Annotation/geneid2go.map")
  for(ensid in unique(biomartResults[,"ensembl_gene_id"])){
    idxes <- which(biomartResults[,"ensembl_gene_id"] == ensid)
    goids <- biomartResults[idxes,"go_id"]
    emptygo <- which(goids=="")
    if(length(emptygo) > 0) goids <- goids[-emptygo]
    if(length(goids) > 0) cat(ensid,"\t", paste(goids, collapse=", "),"\n", file="Annotation/geneid2go.map", append=TRUE, sep="")
  }
}

colfunc       <- colorRampPalette(c("white", "black"))

createGOplot <- function(datasubset, GOtable, x, arrays, colfunc, folder){
  dataonly <- as.matrix(datasubset[,arrays[,"AtlasID"]])
  rownames(dataonly) <- datasubset[,"mgi_symbol"]
  colnames(dataonly) <- apply(arrays[,c("Strain","Tissue")], 1, paste0, collapse="_")
  png(paste0(folder, gsub("GO:","",GOtable[x,"GO.ID"]),"-",GOtable[x,"Term"],".png"))
    heatmap(dataonly, col=colfunc(40))
  dev.off()
}

geneID2GO     <- readMappings(file = "Annotation/geneid2go.map")
topDiffGenes  <- function(x){ return(x < 1e-12) }                                                 # High threshold for tissue analysis

### TODO: After we create the datasubset variable, some probes might NOT show the pattern we're interested in, we need to filter for tissue_p < threshold

### Gene ontology of Hypothalamus
genelist          <- alldata[,"Tissue_P"]                                                         # Create a gene list
names(genelist)   <- alldata[,"ensembl_gene_id"]                                                  # Add the names
genelist[which(alldata[,"HT"] < alldata[,"GF"])] <- 1                                             # Which HT >= GF (so put the ones where HT < GT to P=1)

GOdata            <- new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher      <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
GOHT <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)
for(x in 1:nrow(GOHT)){
  allGenesGOHT  <- genesInTerm(GOdata, whichGO = GOHT[x,"GO.ID"])[[1]]
  genesToGOHT   <- names(which(topDiffGenes(genelist)))
  datasubset    <- alldata[which(alldata[,"ensembl_gene_id"] %in% genesToGOHT[which(genesToGOHT %in% allGenesGOHT)]),]
  write.table(datasubset, file=paste0("GO/HT/",gsub("GO:","",GOHT[x,"GO.ID"]),"-",GOHT[x,"Term"],".txt"), sep="\t", row.names=FALSE)
  createGOplot(datasubset, GOHT, x, arrays, colfunc, "GO/HT/")
}

### Gene ontology of Gonadal fat
genelist          <- alldata[,"Tissue_P"]                                                         # Create a gene list
names(genelist)   <- alldata[,"ensembl_gene_id"]                                                  # Add the names
genelist[which(alldata[,"GF"] < alldata[,"HT"])] <- 1                                             # Which GF >= HT (so put the ones where GT < HT to P=1)

GOdata            <- new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher      <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
GOGF <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)
for(x in 1:nrow(GOGF)){
  allGenesGOHT  <- genesInTerm(GOdata, whichGO = GOGF[x,"GO.ID"])[[1]]
  genesToGOHT   <- names(which(topDiffGenes(genelist)))
  datasubset    <- alldata[which(alldata[,"ensembl_gene_id"] %in% genesToGOHT[which(genesToGOHT %in% allGenesGOHT)]),]
  write.table(datasubset, file=paste0("GO/GF/",gsub("GO:","",GOGF[x,"GO.ID"]),"-",GOGF[x,"Term"],".txt"), sep="\t", row.names=FALSE)
  createGOplot(datasubset, GOGF, x, arrays, colfunc, "GO/GF/")
}

### Gene ontology of strain differences
topDiffGenes <- function(x){ return(x < 1e-4) }                                                   # Lower threshold for strain analysis

genelist          <- alldata[,"Strain_P"]                                                         # Create a gene list
names(genelist)   <- alldata[,"ensembl_gene_id"]                                                  # Add the names

GOdata            <- new("topGOdata", ontology = "BP", allGenes = genelist, geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher      <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
GOstrain <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)

for(x in 1:nrow(GOstrain)){
  allGenesGOHT  <- genesInTerm(GOdata, whichGO = GOstrain[x,"GO.ID"])[[1]]
  genesToGOHT   <- names(which(topDiffGenes(genelist)))
  datasubset    <- alldata[which(alldata[,"ensembl_gene_id"] %in% genesToGOHT[which(genesToGOHT %in% allGenesGOHT)]),]
  write.table(datasubset, file = paste0("GO/Strain/", gsub("GO:","",GOstrain[x,"GO.ID"]), "-", GOstrain[x,"Term"], ".txt"), sep="\t", row.names=FALSE)
  createGOplot(datasubset, GOstrain, x, arrays, colfunc, "GO/Strain/")                            # We might want to split this by tissue to make the effects more clear
}
