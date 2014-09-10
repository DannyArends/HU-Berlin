# geneOntology.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
#

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
RPKM <- read.table("Analysis/BFMI_RPKM_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")

if(!file.exists("GeneOntology/GOannotation.txt")){
  library(biomaRt)                                                                                      # Biomart
  bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                              # Biomart for mouse genes
  biomartResults <- NULL
  for(x in seq(1, length(allgenes), 1000)){                                                             # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(allgenes))                                                            # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    res.biomart <- getBM(c("ensembl_gene_id","go_id"),                                                  # Use biomart to retrieve GO terms
                          filters="ensembl_gene_id", values=allgenes[x:xend], mart=bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
  }
  write.table(biomartResults, file="GeneOntology/GOannotation.txt", sep="\t", row.names=FALSE)
}else{
  biomartResults <- read.table("GeneOntology/GOannotation.txt", sep="\t", header=TRUE)
}

if(!file.exists("GeneOntology/geneid2go.map")){                                                                      # Create the ENSG to GO map only if it doesn't exists
  cat("", file="GeneOntology/geneid2go.map")
  for(ensid in unique(biomartResults[,"ensembl_gene_id"])){
    idxes <- which(biomartResults[,"ensembl_gene_id"] == ensid)
    goids <- biomartResults[idxes,"go_id"]
    emptygo <- which(goids=="")
    if(length(emptygo) > 0) goids <- goids[-emptygo]
    if(length(goids) > 0) cat(ensid,"\t", paste(goids, collapse=", "),"\n", file="GeneOntology/geneid2go.map", append=TRUE, sep="")
  }
}

# Do gene ontology
doGO <- function(allgenes, selected){
  genelist <- rep(0, length(allgenes))                                                                # Create a gene list
  names(genelist) <- allgenes                                                                         # Add the names
  genelist[selected] <- 1                                                                             # Set the switched genes to 1

  library(topGO)

  geneID2GO     <- readMappings(file = "GeneOntology/geneid2go.map")
  GOdata        <- new("topGOdata", ontology = "BP", allGenes = as.factor(genelist), annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes        <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)

  #pdf("GeneOntologyTree.pdf")
    showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
  #dev.off()
  return(allRes)
}

switched <- apply(cbind(RPKM[,"A.D_BFMI860.12xB6N"], RPKM[,"A.D_B6NxBFMI860.12"]),1,function(x){        # Which expressions switch from B6N to BFMI (or vise versa) between the different directions
  if(x[1] == "BFMI" && x[2] == "B6N")  return(1)
  if(x[1] == "B6N"  && x[2] == "BFMI") return(1)
  return(0)
})
sum(switched)

alwaysBFMI <- apply(cbind(RPKM[,"A.D_BFMI860.12xB6N"], RPKM[,"A.D_B6NxBFMI860.12"]),1,function(x){      # Which expressions shows BFMI between the different directions
  if(x[1] == "BFMI" && x[2] == "BFMI") return(1)
  return(0)
})
sum(alwaysBFMI)

alwaysB6N <- apply(cbind(RPKM[,"A.D_BFMI860.12xB6N"], RPKM[,"A.D_B6NxBFMI860.12"]),1,function(x){       # Which expressions shows B6N between the different directions
  if(x[1] == "B6N" && x[2] == "B6N") return(1)
  return(0)
})
sum(alwaysB6N)

doGO(RPKM[, "ensembl_gene_id"], RPKM[which(switched == 1), "ensembl_gene_id"])
doGO(RPKM[, "ensembl_gene_id"], RPKM[which(alwaysBFMI == 1), "ensembl_gene_id"])
doGO(RPKM[, "ensembl_gene_id"], RPKM[which(alwaysB6N == 1), "ensembl_gene_id"])

