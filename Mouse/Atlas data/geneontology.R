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

DEstrain    <- which(as.numeric(alldata[,"Strain_P"]) < DEthreshold)    # Update the strain affected probes, since the ordering has changed
DEtissue    <- which(as.numeric(alldata[,"Tissue_P"]) < DEthreshold)    # Update the tissue affected probes, since the ordering has changed

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

tissueDE        <- alldata[DEtissue,]
geneID2GO         <- readMappings(file = "Annotation/geneid2go.map")

genelist          <- rep(0, length(allgenes))                                                         # Create a gene list
names(genelist)   <- allgenes                                                                         # Add the names
upinHT            <- tissueDE[which(as.numeric(tissueDE[,"HT"]) > as.numeric(tissueDE[,"GF"])),"ensembl_gene_id"]
genelist[upinHT]  <- 1

GOdata            <- new("topGOdata", ontology = "BP", allGenes = as.factor(genelist), annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher      <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)

genelist          <- rep(0, length(allgenes))                                                         # Create a gene list
names(genelist)   <- allgenes                                                                         # Add the names
upinGF            <- tissueDE[which(as.numeric(tissueDE[,"GF"]) > as.numeric(tissueDE[,"HT"])),"ensembl_gene_id"]
genelist[upinGF]  <- 1

GOdata            <- new("topGOdata", ontology = "BP", allGenes = as.factor(genelist), annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher      <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)

