#
# Correlation and differential correlation analysis for reciprocal crosses
#

## Functions
tissueHeatmap <- function(alldata, tissue = "Q", save.results = FALSE) {
  rr <- which(alldata[,paste0("pValues_", tissue)] < 0.05 & (alldata[,paste0("ratios_", tissue)] < 0.83 | alldata[,paste0("ratios_", tissue)] > 1.2))
  cc <- grep(paste0("_", tissue), colnames(alldata))
  cc <- cc[cc < 35]

  correlations <- NULL
  for(x in rr){
    correlations <- rbind(correlations, cor(as.numeric(alldata[x,cc]), t(alldata[rr,cc])))
  }
  rownames(correlations) <- colnames(correlations)
  if(save.results) write.table(correlations, file=paste0("Correlations_", tissue, ".txt"), sep="\t")
  heatmap(correlations, main=paste0("Tissue: ", tissue), scale='none')
}

plotVersus <- function(alldata, gene1 = "ENSMUSG00000065629", gene2 = "ENSMUSG00000104827") {
  expdata <- alldata[,8:37]

  cols <- as.numeric(grepl("^B6N_", colnames(expdata)))
  cols <- cols + (as.numeric(grepl("^BFMI860.12_", colnames(expdata))) * 2)
  cols <- cols + (as.numeric(grepl("^BFMI860.12xB6N", colnames(expdata))) * 3)
  cols <- cols + (as.numeric(grepl("^B6NxBFMI860.12", colnames(expdata))) * 4)
  
  matBFMI <- which(grepl("^BFMI860.12xB6N", colnames(expdata)))
  matB6N <- which(grepl("^B6NxBFMI860.12", colnames(expdata)))
  
  rA <- cor(t(expdata[c(gene1, gene2), matBFMI]))[2,1]
  rB <- cor(t(expdata[c(gene1, gene2), matB6N]))[2,1]  

  r <- cor(t(expdata[c(gene1, gene2), ]))[2,1]

  plot(t(expdata[c(gene1, gene2), ]), col=c("gray", "orange", "blue", "green")[cols], pch=19)
  abline(a = 0, b = r, lty = 2, lwd=3)                          # Estimate intercept:  intercept = mean of Y - Slope * mean of X
  abline(a = 0, b = rA, lty = 2, col="orange", lwd=2)
  abline(a = 0, b = rB, lty = 2, col="gray")
  points(t(expdata[c(gene1, gene2), ]), col=c("gray", "orange", "blue", "green")[cols], pch=19)
}

doCorrelation <- function(expdata, group, output.name, save.results = FALSE, plot.results = FALSE){
  against <- t(expdata[, group])
  outfile <- paste0(output.name,".txt")
  if(!file.exists(outfile)){
    if(save.results) cat(paste(colnames(against), collapse="\t"), "\n", sep="", file = outfile)
    for(x in 1:ncol(against)){
      name <- colnames(against)[x]
      correlations <- round(cor(as.numeric(against[,x]), against, method="spearman"), d = 2)
      if(save.results) cat(paste0(name, "\t", paste0(correlations, collapse="\t")), "\n", sep="", file = outfile, append = TRUE)
    }
  }else{ cat("Loading results from disk\n"); }
  corMatrix <- read.csv(outfile, sep = "\t", header=TRUE, row.names=1)
  if(plot.results) {
    if(save.results) png(paste0(output.name,".png"), width=1024, height=1024)
      op <- par(mai = c(5,5,5,5))
      heatmap(corMatrix, col=c("red", "white", "blue"), breaks=c(-1, -0.5, 0.5, 1), scale = "none", Colv=NA, Rowv=NA, main = output.name, cexRow = 2, cexCol = 2, margins=c(12,8))
    if(save.results) dev.off()
  }
  return(corMatrix)
}

createNetwork <- function(corDIFF, above = 1.5) {
  corD <- abs(corDIFF) > above
  graph <- NULL
  for(x in 1:(nrow(corD)-1)) {
    ii <- which(corD[x,(x+1): ncol(corD)])
    gname <- as.character(expdata[rownames(corD)[x],"gene_name"])
    for(i in ii){
      graph <- rbind(graph, c(rownames(corD)[x], rownames(corD)[i]))
    }
  }
  return(graph)
}

source("D:/Github/HU-Berlin/Mouse/RNASeq2/graph.R")

getEnsembleIDinPathway <- function(patwaydata, name = "Endocytosis") {
  pathways <- as.character(patwaydata[,1])
  mpi <- which(pathways == name)

  ensembleIDs <- unlist(lapply(strsplit(strsplit(as.character(patwaydata[mpi, ncol(patwaydata)]), ";")[[1]],"|",fixed=TRUE),"[",3))
  ensembleIDs <- ensembleIDs[!is.na(ensembleIDs)]
  return(ensembleIDs)
}

plotPathways <- function(patwaydata, corBFMI, corB6N, corDIFF, name1 = "Endocytosis", name2 = "Circadian rhythm", hide.diff = 1, orderByDiff = TRUE) {
  ensID1 <- getEnsembleIDinPathway(patwaydata, name1)                   # Pathway 1
  ensID2 <- getEnsembleIDinPathway(patwaydata, name2)                   # Pathway 2
  cat("Genes in pathways:", length(ensID1), "&", length(ensID2), "\n")
  ensID <- c(ensID1, ensID2)                                            # Combine the pathways

  iD <- corDIFF[ensID, ensID]
  iD <- apply(iD, 1, as.numeric)
  nptIens <- colnames(iD)[which(apply(iD, 1, function(x){all(abs(x) < hide.diff)}))]        # Not interesting (since they do not show a high dCOR)
  if(length(nptIens) > 0){
    ensID1 <- ensID1[-which(ensID1 %in% nptIens)]
    ensID2 <- ensID2[-which(ensID2 %in% nptIens)]
  }
  cat("Genes in pathways (after filtering):", length(ensID1), "&", length(ensID2), "\n")
  
  if(orderByDiff){
    ensID1 <- ensID1[hclust(dist(corDIFF[ensID1,]))$order]
    ensID2 <- ensID2[hclust(dist(corDIFF[ensID2,]))$order]
  }
  
  ensID <- c(ensID1, ensID2)

  iD <- corDIFF[ensID, ensID]
  iD <- apply(iD, 1, as.numeric)
  rownames(iD) <- expdata[colnames(iD),"gene_name"]; colnames(iD) <- expdata[colnames(iD),"gene_name"]

  iBFMI <- corBFMI[ensID, ensID]
  iBFMI <- apply(iBFMI, 1, as.numeric)
  rownames(iBFMI) <- expdata[colnames(iBFMI),"gene_name"]; colnames(iBFMI) <- expdata[colnames(iBFMI),"gene_name"]

  iB6N <- corB6N[ensID, ensID]
  iB6N <- apply(iB6N, 1, as.numeric)
  rownames(iB6N) <- expdata[colnames(iB6N),"gene_name"];colnames(iB6N) <- expdata[colnames(iB6N),"gene_name"]

  colz <- c("tomato4", "tomato3", "tomato2", "tomato1", "white", "slategray1", "slategray2", "slategray3", "slategray4")
  breakz <- c(-1, -0.9, -0.8, -0.6, -0.5, 0.5, 0.6, 0.8, 0.9, 1)
  
  x <- (length(ensID1)-0.5) / (ncol(iBFMI)-1)
  y <- (length(ensID2)-0.5) / (ncol(iBFMI)-1)

  
  op <- par(mfrow=c(1,3), mar=c(10,8,4,2) + 0.1)
  image(iBFMI, col = colz, breaks = breakz, xaxt='n', yaxt='n', main="matBFMI")
  box()
  axis(1, at=seq(0, ncol(iBFMI)-1) / (ncol(iBFMI)-1), rownames(iBFMI), las=2)
  axis(2, at=seq(0, ncol(iBFMI)-1) / (ncol(iBFMI)-1), rownames(iBFMI), las=2)
  abline(v=(length(ensID1)-0.5) / (ncol(iBFMI)-1))
  abline(h=(length(ensID1)-0.5) / (ncol(iBFMI)-1))
  mtext(c(name1,name2),1,5, at = c(x/2, x + 0.5*y),cex=0.7)
  mtext(c(name1,name2),2,5, at = c(x/2, x + 0.5*y),cex=0.7)

  image(iB6N, col = colz, breaks = breakz, xaxt='n', yaxt='n', main="matB6N")
  box()
  axis(1, at=seq(0, ncol(iB6N)-1) / (ncol(iB6N)-1), rownames(iB6N), las=2)
  axis(2, at=seq(0, ncol(iB6N)-1) / (ncol(iB6N)-1), rownames(iB6N), las=2)
  abline(v=(length(ensID1)-0.5) / (ncol(iB6N)-1))
  abline(h=(length(ensID1)-0.5) / (ncol(iB6N)-1))
  mtext(c(name1,name2),1,5, at = c(x/2, x + 0.5*y),cex=0.7)
  mtext(c(name1,name2),2,5, at = c(x/2, x + 0.5*y),cex=0.7)

  image(iD, col = c("tomato3", "white", "slategray2"), breaks = c(-2, -hide.diff, hide.diff, 2), xaxt='n', yaxt='n', main="Difference: matBFMI - matB6N")
  box()
  axis(1, at=seq(0, ncol(iD)-1) / (ncol(iD)-1), rownames(iD), las=2)
  axis(2, at=seq(0, ncol(iD)-1) / (ncol(iD)-1), rownames(iD), las=2)
  abline(v=(length(ensID1)-0.5) / (ncol(iD)-1))
  abline(h=(length(ensID1)-0.5) / (ncol(iD)-1))
  mtext(c(name1,name2),1,5, at = c(x/2, x + 0.5*y),cex=0.7)
  mtext(c(name1,name2),2,5, at = c(x/2, x + 0.5*y),cex=0.7)
}

# Load input RPKM data from RNA-Seq experiment
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross/Analysis")
alldata <- read.table(file=paste0("RPKM_norm_log_stats.txt"), sep="\t", header=TRUE)
rownames(alldata) <- alldata[,"ensembl_gene_id"]

matBFMI <- which(grepl("^BFMI860.12xB6N", colnames(alldata)))
matB6N <- which(grepl("^B6NxBFMI860.12", colnames(alldata)))

## Per tissue correlation ( Might want to use it for the groups )
tissueHeatmap(alldata, "Q")
tissueHeatmap(alldata, "L")
tissueHeatmap(alldata, "G")

## Gene versus Gene plots
plotVersus(alldata, "ENSMUSG00000065629", "ENSMUSG00000104827")
plotVersus(alldata, "ENSMUSG00000065629", "ENSMUSG00000026072")
plotVersus(alldata, "ENSMUSG00000065629", "ENSMUSG00000025981")

# Reduce the dataset to only have genes which have expression in all tissues
good <- which(apply(alldata[,c(matBFMI,matB6N)], 1, function(x){ all(x > 0.5) }))
expdata <- alldata[good,]

# Create a name conversion table
genename <- function(expdata, name = "Per3") { return(rownames(expdata)[which(expdata[,"gene_name"] == name)]); }

# Load the correlations
corBFMI <- doCorrelation(expdata, matBFMI, "matBFMI_correlation")
corB6N  <- doCorrelation(expdata, matB6N, "matB6N_correlation")
corDIFF <- corBFMI - corB6N

# Which genes show high differential correlation
nAbove <- apply(corDIFF, 1,function(x){length(which(abs(x) > 1.5))})
interesting <- as.numeric(nAbove > 1)
names(interesting) <- rownames(corDIFF)

# BIOMART annotation of the MGI_Description
if(!file.exists("E:/Mouse/DNA/Annotation/biomart/BiomartAnnotation.txt")){
  library(biomaRt)
  ensemblIDs <- as.character(dp4datafull[,"gene_id"])
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")                                                # Biomart for mouse genes
  biomartResults <- NULL
  for(x in seq(0, length(ensemblIDs), 1000)){                                                                        # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(ensemblIDs))                                                                       # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")

    res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("ensembl_gene_id"), 
                        values = ensemblIDs[x:xend], mart = bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
    Sys.sleep(1)
  }
  write.table(biomartResults, file="E:/Mouse/DNA/Annotation/biomart/BiomartAnnotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("E:/Mouse/DNA/Annotation/biomart/BiomartAnnotation.txt", sep="\t", header=TRUE, colClasses=c("character"))
}

mm <- NULL
for(x in names(interesting)){
  symbol <- biomartResults[which(biomartResults[,1] == x),"mgi_symbol"]
  if(length(symbol) == 0) symbol = "-"
  if(length(symbol) > 1) symbol = symbol[1]
  mm <- rbind(mm, as.character(c(x,symbol, interesting[x])))
}

write.table(mm, file="scores_diffcor_stefan.txt", row.names=FALSE, quote=FALSE, sep="\t")

connections <- createNetwork(corDIFF)
groups <- traverseGraph(connections)

# Create input for innateDB, we write out groups with more then 8 members
ii <- which(unlist(lapply(groups,length)) >= 8)
for(x in 1:length(ii)){
  cat(groups[[ii[x]]], sep="\n", file=paste0("group", ii[x], ".txt"))
}

# Load the pathway ORA results from innateDB
pathwaydataGroup1 <- read.csv("group1pathwayORA.txt", sep="\t", header=TRUE, colClasses="character")

# Create some plots of different pathways
plotPathways(pathwaydataGroup1, corBFMI, corB6N, corDIFF, "Endocytosis", "Circadian rhythm")
plotPathways(pathwaydataGroup1, corBFMI, corB6N, corDIFF, "Metabolism", "Circadian rhythm")
plotPathways(pathwaydataGroup1, corBFMI, corB6N, corDIFF, "Metabolism", "Endocytosis")
plotPathways(pathwaydataGroup1, corBFMI, corB6N, corDIFF, "Endocytosis", "Oxidative phosphorylation")
plotPathways(pathwaydataGroup1, corBFMI, corB6N, corDIFF, "Circadian rhythm", "Oxidative phosphorylation")

checkData <- function(x, y){
  cat(corBFMI[rownames(expdata)[which(expdata[,"gene_name"] == x)], rownames(expdata)[which(expdata[,"gene_name"] == y)]], "\n")
  cat(corB6N[rownames(expdata)[which(expdata[,"gene_name"] == x)], rownames(expdata)[which(expdata[,"gene_name"] == y)]], "\n")
  cat(corDIFF[rownames(expdata)[which(expdata[,"gene_name"] == x)], rownames(expdata)[which(expdata[,"gene_name"] == y)]], "\n")
}

checkData("Creb1", "Per3")
checkData("Ccnd3", "Anapc7")
checkData("Per3", "Per1")
checkData("Per3", "Polr2h")
