# geneOntology.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
#

library(biomaRt)                                                                                        # Biomart package
library(topGO)                                                                                          # topGO package
  
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
RPKM <- read.table("Analysis/BFMI_RPKM_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character", check.names=FALSE)

allgenes <- RPKM[, "ensembl_gene_id"]

if(!file.exists("GeneOntology/GOannotation.txt")){
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
  cat("Loading biomaRt gene ontology annotation from disk\n")
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

  geneID2GO     <- readMappings(file = "GeneOntology/geneid2go.map")
  GOdata        <- new("topGOdata", ontology = "BP", allGenes = as.factor(genelist), annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes        <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)

  #pdf("GeneOntologyTree.pdf")
    showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
  #dev.off()
  return(allRes)
}

switched <- apply(cbind(RPKM[,"A/D_BFMI860-12xB6N"], RPKM[,"A/D_B6NxBFMI860-12"]),1,function(x){        # Which expressions switch from B6N to BFMI (or vise versa) between the different directions
  if(x[1] == "BFMI" && x[2] == "B6N")  return(-1)
  if(x[1] == "B6N"  && x[2] == "BFMI") return(1)
  return(0)
})
sum(switched)

alwaysBFMI <- apply(cbind(RPKM[,"A/D_BFMI860-12xB6N"], RPKM[,"A/D_B6NxBFMI860-12"]),1,function(x){      # Which expressions shows BFMI between the different directions
  if(x[1] == "BFMI" && x[2] == "BFMI") return(1)
  return(0)
})
sum(alwaysBFMI)

alwaysB6N <- apply(cbind(RPKM[,"A/D_BFMI860-12xB6N"], RPKM[,"A/D_B6NxBFMI860-12"]),1,function(x){       # Which expressions shows B6N between the different directions
  if(x[1] == "B6N" && x[2] == "B6N") return(1)
  return(0)
})
sum(alwaysB6N)

LD1 <- RPKM[,c("F1-V-1004_L", "F1-V-1016_L", "F1-V-1020_L")]                                                        # BFMI cross BFMI860-12xB6N (D1)
LD2 <- RPKM[,c("F1-V-1000_L", "F1-V-1008_L", "F1-V-1012_L")]                                                        # BFMI cross B6NxBFMI860-12 (D2)
pval <- apply(cbind(LD1,LD2),1,function(x){ return(t.test(as.numeric(x[1:3]), as.numeric(x[4:6]))$p.value)})        # T-test for differences

goSwitched <- doGO(RPKM[, "ensembl_gene_id"], RPKM[which(switched != 0), "ensembl_gene_id"])
write.table(goSwitched, "GO_Switched.txt", sep="\t", row.names=FALSE, quote=FALSE)
goMaternal <- doGO(RPKM[, "ensembl_gene_id"], RPKM[which(switched == -1), "ensembl_gene_id"])
write.table(goMaternal, "GO_Maternal.txt", sep="\t", row.names=FALSE, quote=FALSE)
goPaternal <- doGO(RPKM[, "ensembl_gene_id"], RPKM[which(switched == 1), "ensembl_gene_id"])
write.table(goPaternal, "GO_Paternal.txt", sep="\t", row.names=FALSE, quote=FALSE)
goBFMI     <- doGO(RPKM[, "ensembl_gene_id"], RPKM[which(alwaysBFMI == 1), "ensembl_gene_id"])
write.table(goBFMI, "GO_BFMI.txt", sep="\t", row.names=FALSE, quote=FALSE)
goB6N      <- doGO(RPKM[, "ensembl_gene_id"], RPKM[which(alwaysB6N == 1), "ensembl_gene_id"])
write.table(goB6N, "GO_B6N.txt", sep="\t", row.names=FALSE, quote=FALSE)
goDE       <- doGO(RPKM[, "ensembl_gene_id"], RPKM[which(pval < 0.005), "ensembl_gene_id"])
write.table(goDE, "GO_DiffExp_0.005.txt", sep="\t", row.names=FALSE, quote=FALSE)

write.table(RPKM[which(switched == 1), ],   "Expression_Switched.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(RPKM[which(alwaysBFMI == 1), ], "Expression_BFMI.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(RPKM[which(alwaysB6N == 1), ],  "Expression_B6N.txt", sep="\t", row.names=FALSE, quote=FALSE)
write.table(RPKM[which(pval < 0.005), ],    "Expression_Differential_0.005.txt", sep="\t", row.names=FALSE, quote=FALSE)



allGenesGO <- as.character(biomartResults[which(biomartResults[,"go_id"] == goSwitched[1,"GO.ID"]),"ensembl_gene_id"])
allGenesGO

RPKM[which(switched != 0),][which(as.character(RPKM[which(switched != 0), "ensembl_gene_id"] %in% allGenesGO),]

