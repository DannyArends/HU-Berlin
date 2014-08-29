# Functions for the script for Futterungs Versuch
# (c) Nov 2013 Danny Arends
# Make sure express[.exe] is on your path

# Plot number 2 from the DESeq2 package (Unused)
plotF2 <-function(res){
  plot(res$baseMean, pmin(-log10(res$pvalue), 50), log="x", xlab="mean of normalized counts", ylab=expression(-log[10](pvalue)))
  abline(v=10,col="red",lwd=1)
}

# Plot number 4 from the DESeq2 package (Unused)
plotF4 <-function(dds){
  rld <- rlogTransformation(dds, blind=FALSE)
  par(mfrow=c(1,2))
  notAllZero <- (rowSums(counts(dds))>0)
  meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim = c(0,2.5))
  meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))
}
# Function filter the results and adjust the pvalues after filtering (Unused)
filterResults <- function(res, significance = 0.05, cutoff = 10, adjustPvalues = TRUE){
  use <- res$baseMean >= cutoff & !is.na(res$pvalue)                          # Filter the output results
  resFilt <- res[use,]
  if(adjustPvalues) resFilt$padj <- p.adjust(resFilt$pvalue, method="BH")     # Adjust the P-values after filtering
  cat("Number of significant (Before filtering): ", sum(res$padj < significance, na.rm=TRUE),"\n")
  cat("Number of significant (After filtering): ", sum(resFilt$padj < significance, na.rm=TRUE),"\n")
  return(resFilt)
}

# Function to visualize the multiple testing adjustment (Unused)
plotMultiTestAdjustment <- function(results, main, show = 0.01, alpha = 0.1){
  orderInPlot <- order(results$pvalue)
  showInPlot <- (results$pvalue[orderInPlot] <= show)
  plot(seq(along=which(showInPlot)), results$pvalue[orderInPlot][showInPlot], pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]), main=main)
  abline(a=0, b=alpha/length(results$pvalue), col="red3", lwd=2)
}

# Function that removes all the non significant pvalues
# It also transforms the Log2FoldChange into FoldChange for clarity
onlySignificant <- function(model, significance = 0.05){ 
  subSet <- as.matrix(model[which(model[,"pvalue"] <= significance), -c(3,4)])
  subSet[,"log2FoldChange"] <- sign(subSet[,"log2FoldChange"]) * (2^abs(subSet[,"log2FoldChange"]))
  colnames(subSet) <- c("baseMean", "FoldChange","pvalue","padj")
  return(subSet)
}

# ------------------------ Main analysis function around DESeq2 ------------------------ #
# return a LIST of [[2]] elements
# [[1]] = The model
# [[2]] = The normalized counts
doDifferentialExpression <- function(countMatrix, colData, design, orderPvals = TRUE, na.rm = TRUE){
  dds <- DESeqDataSetFromMatrix(countData = countMatrix, colData = as.data.frame(colData), design = design)
  colData(dds)$tissue <- factor(colData(dds)$tissue, levels=c("Jejunum", "Ileum"))
  colData(dds)$celltype <- factor(colData(dds)$celltype, levels=c("PP", "LN"))
  colData(dds)$condition <- factor(colData(dds)$condition, levels=c("untreated","treated"))

  dds <- DESeq(dds)
  counts <- counts(dds,normalized=TRUE)
  res <- results(dds)

  if(orderPvals) res <- res[order(res$padj),]
  if(na.rm) res <- res[!is.na(res$pvalue),]
  return(list(res, counts))
}

# Function to get the normalized counts out of the model
getNormalizedCounts <- function(model){ return(round(model[[2]][rownames(onlySignificant(model[[1]] )),], d = 1)) }

# Function to download data from biomart and create an annotation matrix
createAnnotationMatrix <- function(genes, biomart = "ensembl", dataset="sscrofa_gene_ensembl"){
  cat("Start requesting data from biomart: ", biomart, "-", dataset,"\n")
  mart <- useMart(biomart = biomart, dataset = dataset)
  
  resultsHGNC <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","name_1006","description"), filters = "hgnc_symbol", values = genes, mart = mart)
  resultsEGID <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","name_1006","description"), filters = "ensembl_gene_id", values = genes, mart = mart)
  results <- rbind(resultsEGID,resultsHGNC)
  cat("Retrieved: ", nrow(results), "records\n")
  annotations <- lapply(unique(results[,1]),function(x){
    ids <- which(results[,1] == x)
    c(x, results[ids[1],2], paste(unlist(results[ids,"name_1006"]), collapse="; "), results[ids[1], "description"])
  })

  cat("Building AnnotationMatrix\n")
  annotationM <- NULL
  for(x in 1:length(annotations)){
    annotationM <- rbind(annotationM, c(annotations[[x]][1], annotations[[x]][2], annotations[[x]][3], annotations[[x]][4]))
    if(x %% 1000 == 0) cat("Done",x,"/",length(annotations),"\n")
  }
  invisible(annotationM)
}

# Annotate a model using the annotationMatrix
# We query twice (Which is not needed any more because we now have good identifiers)
annotateModel <- function(model, annotationM){
  identifiers <- unlist(lapply(strsplit(rownames(model), ":"),"[", 1))
  annotatedModel <- NULL
  for(x in 1:length(identifiers)){
    id <- identifiers[x]
    annotation  <- "NA"
    description <- "NA"
    if((id %in% annotationM[,1])){
      annotation  <- annotationM[which(annotationM[,1] == id),3]
      description <- annotationM[which(annotationM[,1] == id),4]
    }else if((id %in% annotationM[,2])){
      annotation  <- annotationM[which(annotationM[,2] == id),3]
      description <- annotationM[which(annotationM[,2] == id),4]
    }else{
      #cat("[Warning] ", id, "is not found in the annotation\n")
    }
    if(annotation == "") annotation  <- "NA"
    annotatedModel <- rbind(annotatedModel, c(model[x,], annotation, description))
  }
  colnames(annotatedModel) <- c(colnames(model),"GO Annotation","Description")
  rownames(annotatedModel) <- rownames(model)
  additionalMatrix <- strsplit(rownames(annotatedModel), ":")
    geneNames <- unlist(lapply(additionalMatrix,"[", 2))
    geneID <- unlist(lapply(additionalMatrix,"[", 1))
    transcriptID <- unlist(lapply(additionalMatrix,"[", 3))
    firstColumns <- cbind(geneNames, geneID, transcriptID)
  annotatedModel <- as.data.frame(annotatedModel)
  annotatedModel[,1:length(colnames(model))] <- apply(annotatedModel[,1:length(colnames(model))],2,as.numeric) # Force numeric to prevent excel from messing up
  annotatedModel <- cbind(firstColumns, annotatedModel)
  return(annotatedModel)
}

# Get the names of the genes which have a up / down regulation above a certain fold change
getUpRegulated   <- function(model, cutoff = 1.25){ names(which(model[,"FoldChange"] >= cutoff)) }
getDownRegulated <- function(model, cutoff = 1.25){ names(which(model[,"FoldChange"] <= cutoff)) }

# Create groups from the model data which are either Up or Down Regulated (Use the above functions as input to FUN)
createGroup <- function(modelData, FUN = getUpRegulated){ return(unlist(lapply(strsplit(FUN(modelData), "_"),"[", 1))) }

# Create the Venn Diagram using a list of models, where we apply FUN to extract the Up / Down regulated
createVennDiagram <- function(ListOfModels, FUN = getUpRegulated){
  groups <- vector("list", length(ListOfModels))
  modelID <- 1
  for(model in ListOfModels){
    groups[[modelID]] <- unique(createGroup(model, FUN))
    modelID <- modelID + 1
  }
  venn(groups)
  invisible(groups)
}

# Create a nice excel sheet with the genes that are UP/DOWN regulated
# Assumes the length of ListOfGroups == 7 (Because it adds column headers to it)
createExcelSheet <- function(ListOfGroups){
  cmax <- 0
  for(group in ListOfGroups){ cmax <- max(cmax, length(group)) }    # Figure out which list is the longest so we can create a matrix
  
  toExcel <- matrix("", cmax, length(ListOfGroups))
  modelID <- 1
  for(group in ListOfGroups){
    toExcel[1:length(group), modelID] <- unlist(lapply(strsplit(group, ":"),"[", 2))
    reNameToENSSSCG <- which(toExcel[1:length(group), modelID] == "NA")                                 # Gene names can be NA
    toExcel[reNameToENSSSCG, modelID] <- unlist(lapply(strsplit(group[reNameToENSSSCG], ":"),"[", 1))   # Use the ENSSSCG in that case
    modelID <- modelID + 1
  }
  colnames(toExcel) <- c("Overlap All", "Overlap PP & LN", "in PP not in LN", "in LN not in PP", "Overlap JJ & IL", "in JJ not in IL", "in IL not in JJ")
  return(toExcel)
}
