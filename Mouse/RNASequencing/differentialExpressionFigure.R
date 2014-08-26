# differentialExpressionFigure.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
#
# Create a figure for the RNA sequencing data (pre-processed by MDC)

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
  
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/MPI_RPKM_ANALYSIS/")
RPKM <- read.table("BFMI_RPKM_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")

mlength <- max(chrInfo[,"Length"])

#png("MaternalOrigin.png", width = 2000, height = 1000)
  #op <- par(cex = 2.5)
  plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="Dominant maternal origin", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")

  aa <- apply(RPKM, 1,function(x){
    yloc <- match(as.character(x["chromosome_name"]), chromosomes); xloc <- as.numeric(x["start_position"])
    #ttest <- as.numeric(x["tTest.1"])
    #if(!is.na(ttest) && ttest < 0.1){
      col <- "white"
      if(x["A.D_BFMI860.12xB6N"] == "B6N") col <- "gray"
      if(x["A.D_BFMI860.12xB6N"] == "BFMI") col <- "orange"
      if(x["A.D_BFMI860.12xB6N"] == "ADDITIVE") col <- "white"
      if(col != "white") points(x=xloc, y=yloc + 0.1, pch=15, col=col,cex=0.9)
      col <- "white"
      if(x["A.D_B6NxBFMI860.12"] == "B6N") col <- "gray"
      if(x["A.D_B6NxBFMI860.12"] == "BFMI") col <- "orange"
      if(x["A.D_B6NxBFMI860.12"] == "ADDITIVE") col <- "white"
      if(col != "white") points(x=xloc, y=yloc - 0.1, pch=15, col=col,cex=0.9)
    #}
  })

  cnt <- 1
  aa <- apply(chrInfo,1,function(x){
    lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1,lwd=2)
    cnt <<- cnt + 1
  })

  axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
  axis(1, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7)
  legend("topright", c("BFMI like", "B6N like"), fill=c("orange","gray"))
#dev.off()

cat("Maternal BFMI: B6N:", sum(RPKM[,"A.D_BFMI860.12xB6N"] == "B6N"), "BFMI:", sum(RPKM[,"A.D_BFMI860.12xB6N"] == "BFMI"),"\n")
cat("Maternal B6: B6N:", sum(RPKM[,"A.D_B6NxBFMI860.12"] == "B6N"), "BFMI:", sum(RPKM[,"A.D_B6NxBFMI860.12"] == "BFMI"),"\n")

# Continue with the Maternal BFMI
which(RPKM[,"A.D_B6NxBFMI860.12"] == "B6N")

selected <- RPKM[which(RPKM[,"A.D_BFMI860.12xB6N"] == "BFMI"), "ensembl_gene_id"]
allgenes <- RPKM[, "ensembl_gene_id"]

geneList <- rep(0, length(allgenes))
names(geneList) <- allgenes
geneList[selected] <- 1

if(!file.exists("GOannotation.txt")){
  library(biomaRt)                                                                                                    # Biomart
  bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
  biomartResults <- NULL
  for(x in seq(1, length(allgenes), 1000)){                                                             # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(allgenes))                                                            # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    res.biomart <- getBM(c("ensembl_gene_id","go_id"),                                                  # Use biomart to retrieve GO terms
                          filters="ensembl_gene_id", values=allgenes[x:xend], mart=bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
  }
  write.table(biomartResults, file="GOannotation.txt", sep="\t", row.names=FALSE)
}else{
  biomartResults <- read.table("GOannotation.txt", sep="\t", header=TRUE)
}

cat("", file="geneid2go.map")
for(ensid in unique(biomartResults[,"ensembl_gene_id"])){
  idxes <- which(biomartResults[,"ensembl_gene_id"] == ensid)
  goids <- biomartResults[idxes,"go_id"]
  emptygo <- which(goids=="")
  if(length(emptygo) > 0) goids <- goids[-emptygo]
  if(length(goids) > 0) cat(ensid,"\t", paste(goids, collapse=", "),"\n", file="geneid2go.map", append=TRUE, sep="")
}

library(topGO)

geneID2GO     <- readMappings(file = "geneid2go.map")
GOdata        <- new("topGOdata", ontology = "BP", allGenes = as.factor(geneList), annot = annFUN.gene2GO, gene2GO = geneID2GO)
resultFisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
allRes        <- GenTable(GOdata, classicFisher = resultFisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)

pdf("GeneOntologyTree.pdf")
  showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
dev.off()

