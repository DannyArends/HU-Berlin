# geneOntology.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Oct, 2014
# first written Oct, 2014
#

library(biomaRt)                                                                                        # Biomart package
library(topGO)                                                                                          # topGO package
  
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
RPKM <- read.table("Analysis/BFMI_RPKM_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character", check.names=FALSE)

allgenes        <- RPKM[, "ensembl_gene_id"]
genelist        <- rep(0, length(allgenes))                                                      # Create a gene list
names(genelist) <- allgenes                                                                      # Add the names
genelist[1]     <- 1                                                                             # Set a random gene to 1, (the first one)
 
geneID2GO     <- readMappings(file = "GeneOntology/geneid2go.map")
GOdata        <- new("topGOdata", ontology = "BP", allGenes = as.factor(genelist), annot = annFUN.gene2GO, gene2GO = geneID2GO)

GetAllTermsInGO <- function(RPKM, GOdata, whichGO = ""){
  allGenesGO    <- genesInTerm(GOdata, whichGO = whichGO)[[1]]
  return(RPKM[which(as.character(RPKM[, "ensembl_gene_id"]) %in% allGenesGO),])
}

GO0006631 <- GetAllTermsInGO(RPKM, GOdata, "GO:0006631")