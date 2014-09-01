# Script for Futterungs Versuch
# (c) Nov 2013 Danny Arends
# NOTE: Make sure express[.exe] is on your path

# ------------------------ Set up the environment ------------------------ #

setwd("D:/Github/Gudrun")        # Location of script files
source("functions.R")     # Source the functions

# If you need to install the packages please uncomment the following lines
#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("vsn")
#biocLite("biomaRt")
#install.packages("xlsx")
#install.packages("gplots")
#install.packages("extrafont")     ## Extra packages to be able to use Arial in plots

library("extrafont")
library("parallel")
library("DESeq2")
library("vsn")
library("xlsx")
library("biomaRt")
library("gplots")

font_import()

setwd("H:/Expressionsanalyse/results_pig_Ensembl73") # Location of input files

# ------------------------ Names of the samples provided ------------------------ #
sampleNames <- c("JE_PP_c", "JE_LN_c", "Il_PP_c","Il_LN_c","JE_PP_t","JE_LN_t","Il_PP_t","Il_LN_t")
nSamples    <- length(sampleNames)

# ------------------------ Express alignment to reference genome ------------------------ #
# Once aligned we don't need to re-do the alignment, redo the alignment by commenting out these lines
#
#cl <- makeCluster(getOption("cl.cores", 4))
#parLapply(cl, 1:nSamples, function(x){
#  reference <- "Sus_scrofa_10.2_v73_transcripts.fa"
#  logfile   <- paste("log", x, ".txt", sep="")
#  cmd       <- paste("express --logtostderr -o Alignment",x," ",reference," DNA10447-",x,"-L-4_preprocessed.srt.bam 2> ",logfile, sep="")
#  shell(cmd)
#})
#stopCluster(cl)

# ------------------------ Load RAW Count Matrix ------------------------ #
baseMatrix <- read.csv("Alignment1/results.xprs", sep="\t", row.names = 2) # Used as row reference, we use the ordering of rows in file 1
countMatrix <- NULL
for(x in 1:nSamples){
  alignmentResults <- read.csv(paste("Alignment",x,"/results.xprs",sep=""), sep="\t", row.names=2)
  # (Re-)order the loaded alignmentResults to baseMatrix
  countMatrix <- cbind(countMatrix, alignmentResults[match(rownames(baseMatrix), rownames(alignmentResults)), "tot_counts"])
  cat("Done ", x, "/", nSamples, "\n",sep="")
}
rownames(countMatrix) <- rownames(baseMatrix)
colnames(countMatrix) <- sampleNames

# NOTE: Fix/HACK for DEseq2 package 
# We throw away all the genes which have very low counts 10 in 8 samples = 1.25 read/sample
countMatrix  <- countMatrix[-which(apply(countMatrix,1,sum) < 10),]

# ------------------------ Covariates ------------------------ #
tissue      <- c("Jejunum", "Jejunum", "Ileum", "Ileum", "Jejunum", "Jejunum", "Ileum", "Ileum")
celltype    <- c("PP", "LN", "PP", "LN", "PP", "LN", "PP", "LN")
conditions  <- c("untreated", "untreated", "untreated", "untreated", "treated", "treated", "treated", "treated")
tech        <- rep("single-read", 8)

# ------------------------ Create the column description data for DESeq2 ------------------------ #
colData <- cbind(tissue, celltype, conditions, tech)
colnames(colData) <- c("tissue", "celltype", "condition", "type")
rownames(colData) <- sampleNames

# ------------------------ Download the biomart annotation and create the annotation matrix ------------------------ #
geneNames <- unlist(lapply(strsplit(rownames(countMatrix), ":"),"[", 1))
annotationM <- createAnnotationMatrix(geneNames)

# ------------------------ Model fitting ------------------------ #

# The full model, should perform the best
# 4 replicated per treatment group so P values are semi trustworthy
fullModel    <- doDifferentialExpression(countMatrix, colData, design = ~ tissue + celltype + condition )

# 2 Models for the individual tissues
# Only 2 replicates per treatment group so P values are not to be trusted 
jejunumModel <- doDifferentialExpression(countMatrix[,c(1,2,5,6)], colData[c(1,2,5,6),], design = ~ celltype + condition )
ileumModel   <- doDifferentialExpression(countMatrix[,c(3,4,7,8)], colData[c(3,4,7,8),], design = ~ celltype + condition )

# 2 Models for individual cell types
# Only 2 replicates per treatment group so P values are not to be trusted 
ppModel   <- doDifferentialExpression(countMatrix[,c(1,3,5,7)], colData[c(1,3,5,7),], design = ~ tissue + condition )
lnModel   <- doDifferentialExpression(countMatrix[,c(2,4,6,8)], colData[c(2,4,6,8),], design = ~ tissue + condition )

# Add the normalized and raw counts to the fullModelData
fullModelData <- cbind(getNormalizedCounts(fullModel), onlySignificant(fullModel[[1]]), round(countMatrix[rownames(onlySignificant(fullModel[[1]] )),], d = 1))
colnames(fullModelData)[1:8]   <- sampleNames
colnames(fullModelData)[13:20] <- sampleNames

# Add the normalized and raw counts to the jejunumModelData
jejunumModelData <- cbind(getNormalizedCounts(jejunumModel), onlySignificant(jejunumModel[[1]]), round(countMatrix[rownames(onlySignificant(jejunumModel[[1]] )),c(1,2,5,6)], d = 1))
colnames(jejunumModelData)[1:4]  <- sampleNames[c(1,2,5,6)]
colnames(jejunumModelData)[9:12] <- sampleNames[c(1,2,5,6)]

# Add the normalized and raw counts to the ileumModelData
ileumModelData <- cbind(getNormalizedCounts(ileumModel), onlySignificant(ileumModel[[1]]), round(countMatrix[rownames(onlySignificant(ileumModel[[1]] )),c(3,4,7,8)], d = 1))
colnames(ileumModelData)[1:4]  <- sampleNames[c(3,4,7,8)]
colnames(ileumModelData)[9:12] <- sampleNames[c(3,4,7,8)]

# Add the normalized and raw counts to the ileumModelData
ppModelData <- cbind(getNormalizedCounts(ppModel), onlySignificant(ppModel[[1]]), round(countMatrix[rownames(onlySignificant(ppModel[[1]] )), c(1,3,5,7)], d = 1))
colnames(ppModelData)[1:4]   <- sampleNames[c(1,3,5,7)]
colnames(ppModelData)[9:12]  <- sampleNames[c(1,3,5,7)]

# Add the normalized and raw counts to the ileumModelData
lnModelData <- cbind(getNormalizedCounts(lnModel), onlySignificant(lnModel[[1]]), round(countMatrix[rownames(onlySignificant(lnModel[[1]] )), c(2,4,6,8)], d = 1))
colnames(lnModelData)[1:4]   <- sampleNames[c(2,4,6,8)]
colnames(lnModelData)[9:12]  <- sampleNames[c(2,4,6,8)]

# ------------------------ Add the model data to an excel file ------------------------ #
wb <- createWorkbook()

fullModelSheet <- createSheet(wb, sheetName="fullModel")
addDataFrame(annotateModel(fullModelData, annotationM), fullModelSheet)

jejunumModelSheet <- createSheet(wb, sheetName="jejunumModel")
addDataFrame(annotateModel(jejunumModelData, annotationM), jejunumModelSheet)

ileumModelSheet <- createSheet(wb, sheetName="ileumModel")
addDataFrame(annotateModel(ileumModelData, annotationM), ileumModelSheet)

ppModelSheet <- createSheet(wb, sheetName="ppModel")
addDataFrame(annotateModel(ppModelData, annotationM), ppModelSheet)

lnModelSheet <- createSheet(wb, sheetName="lnModel")
addDataFrame(annotateModel(lnModelData, annotationM), lnModelSheet)

# ------------------------ Create a Venn Diagram for UP regulation ------------------------ #
postscript(file="VennUpRegulation_1.25.eps", width=10, height=10)
  groups <- createVennDiagram(list(fullModelData, jejunumModelData, ileumModelData, ppModelData, lnModelData))
  legend("topleft",c("Fold Change > 1.25", "A = Full Model", "B = Jejunum Model", "C = Ileum Model", "D = Peyer's Patch Model", "E = Lymph Node Model"))
dev.off()

# ------------------------ Create the 'universe' to query the VENN overlaps ------------------------ #
universe <- unique(c(unlist(groups)))
group1.l <- universe %in% groups[[1]]
group2.l <- universe %in% groups[[2]]
group3.l <- universe %in% groups[[3]]
group4.l <- universe %in% groups[[4]]
group5.l <- universe %in% groups[[5]]

# ------------------------ Create the requested groups ------------------------ #

sharedAllUp   <- universe[group1.l & group2.l & group3.l & group4.l& group5.l]  # Shared between all models

sharedPP_LN <- universe[group4.l &  group5.l]                                   # Shared between PP and LN
inPPnotLN   <- universe[group4.l & !group5.l]                                   # In PP but not in LN
inLNnotPP   <- universe[group5.l & !group4.l]                                   # In LN but not in PP

sharedJJ_IL <- universe[group2.l &  group3.l]                                   # Shared between Jejunum and Ileum
inJJnotIL   <- universe[group2.l & !group3.l]                                   # In Jejunum but not in Ileum
inILnotJJ   <- universe[group3.l & !group2.l]                                   # In Ileum but not in Jejunum

# ------------------------ Save the UP regulated groups to Excel  ------------------------ #

toExcel <- createExcelSheet(list(sharedAllUp, sharedPP_LN, inPPnotLN, inLNnotPP, sharedJJ_IL, inJJnotIL, inILnotJJ))

upSheet <- createSheet(wb, sheetName="Up Regulated Genes")
addDataFrame(toExcel, upSheet, row.names=FALSE)

# ------------------------ Create a Venn Diagram for DOWN regulation ------------------------ #
postscript(file="VennDownRegulation_1.25.eps", width=10, height=10)
  groups <- createVennDiagram(list(fullModelData, jejunumModelData, ileumModelData, ppModelData, lnModelData), getDownRegulated)
  legend("topleft",c("Fold Change < -1.25", "A = Full Model", "B = Jejunum Model", "C = Ileum Model", "D = Peyer's Patch Model", "E = Lymph Node Model"))
dev.off()

# ------------------------ Create the universe to query ------------------------ #
universe <- unique(c(unlist(groups)))
group1.l <- universe %in% groups[[1]]
group2.l <- universe %in% groups[[2]]
group3.l <- universe %in% groups[[3]]
group4.l <- universe %in% groups[[4]]
group5.l <- universe %in% groups[[5]]

# ------------------------ Create the requested groups ------------------------ #
sharedAllDown   <- universe[group1.l & group2.l & group3.l & group4.l& group5.l]  # Shared between all models

sharedPP_LN <- universe[group4.l &  group5.l]                                     # Shared between PP and LN
inPPnotLN   <- universe[group4.l & !group5.l]                                     # In PP but not in LN
inLNnotPP   <- universe[group5.l & !group4.l]                                     # In LN but not in PP

sharedJJ_IL <- universe[group2.l & group3.l]                                      # Shared between Jejunum and Ileum
inJJnotIL   <- universe[group2.l & !group3.l]                                     # In Jejunum but not in Ileum
inILnotJJ   <- universe[group3.l & !group2.l]                                     # In Ileum but not in Jejunum

# ------------------------ Save the DOWN regulated groups to Excel  ------------------------ #
toExcel <- createExcelSheet(list(sharedAllDown, sharedPP_LN, inPPnotLN, inLNnotPP, sharedJJ_IL, inJJnotIL, inILnotJJ))

downSheet <- createSheet(wb, sheetName="Down Regulated Genes")
addDataFrame(toExcel, downSheet, row.names=FALSE)

# ------------------------ Save the excel file  ------------------------ #
saveWorkbook(wb, file="FutterungsVersuch.xlsx")

# ------------------------ Create the tissue heatmaps (Up and Down)  ------------------------ #
counts <- getNormalizedCounts(fullModel)
colnames(counts)   <- sampleNames
ratios <- (counts[,5:8] - counts[,1:4])     #Should actually be called DIFFERENCE

scRatios <- t(apply(ratios,1,function(x){ return(rank(x/max(abs(x)))) }))                               # Scale and rank, so we have nice plots

colnames(scRatios) <- c("JE PP","JE LN","IL PP","IL LN")                                                # Use the shorter names
fullNames <- rownames(scRatios)
rownames(scRatios) <- unlist(lapply(strsplit(fullNames, ":"),"[", 2))                                   # Remove chromosome & Location
reNameToENSSSCG <- which(rownames(scRatios) == "NA")                                                    # NA gene names
rownames(scRatios)[reNameToENSSSCG] <- unlist(lapply(strsplit(fullNames[reNameToENSSSCG], ":"),"[", 1)) # Use the ENSSSCG ID in plots

up   <- unlist(lapply(sharedAllUp,    function(x){which(grepl(x, rownames(ratios)))}))                  # Get all the shared UP regulated genes for plotting
down <- unlist(lapply(sharedAllDown,  function(x){which(grepl(x, rownames(ratios)))}))                  # Get all the shared DOWN regulated genes for plotting

# ------------------------ Write the heatmaps to a PDF file  ------------------------ #

postscript(file="UpRegulated_4_tissues.eps", width=10, height=10)
  heatmap(scRatios[c(up),], margins=c(6, 9), main="Up regulated genes", col=gray.colors(4))
  legend("topleft",title="FoldChange", c("High","Medium","Low","Lowest"), col=c(gray.colors(4)),lwd=c(10))
dev.off()

postscript(file="UpRegulated_4_tissues_NOTREE.eps", width=10, height=15)
  op <- par(oma=c(2,2,10,2))
  heatmap(scRatios[c(up),], Colv=NA, margins=c(6, 9), main="Up regulated genes", col=gray.colors(4))
  op <- par(oma=c(0,0,0,0))  
  legend("topleft",title="FoldChange", c("High","Medium","Low","Lowest"), col=c(gray.colors(4)),lwd=c(10))
dev.off()


genenames <- c("ENSSSCG00000015781", "ENSSSCG00000014727", "HBB - ENSSSCT00000036536", "CLU - ENSSSCT00000036649", "HBB - ENSSSCT00000016076", "MZB1", "LAG3", "PTPRCAP", "IGKV-3", "IGKV-5",
  "CCL19", "ENSSSCG00000029057", "ENSSSCG00000020750", "IGKC", "IGKV-11", "IGKJ2", "IGKV-7", "IGKV-6", "ENSSSCG00000008205", "IGLV-5", "CLU - ENSSSCT00000010601", "ENSSSCG00000030271",
  "HEXDC", "IGJ - ENSSSCT00000033582", "IGJ - ENSSSCT00000009792", "IGLC - ENSSSCT00000010999", "IGLC - ENSSSCT00000011005", "ENSSSCG00000010077", "IGLV-8", "IGLV-4", "IGLV-7")

loadfonts(device = "postscript") ## for postscript()

dendro <- heatmap(scRatios[c(down),], margins=c(6, 9), main="Down regulated genes", col=gray.colors(4),Rowv=NA, keep.dendro = TRUE)

postscript(file="DownRegulated_4_tissues.jenny.eps", width=9, height=10, family = "Arial", paper = "special", horizontal=FALSE)
  op <- par(omi = c(0,0,0,2))
  heatmap(scRatios[c(down),][dendro$rowInd,], margins=c(6, 9), main="Down regulated genes", col=gray.colors(4),Colv=NA,Rowv=NA)
  legend("topleft",title="FoldChange", c("High","Medium","Low","Lowest"), col=c(gray.colors(4)),lwd=c(10))
dev.off()

postscript(file="DownRegulated_4_tissues_NOTREE.eps", width=10, height=15)
  op <- par(oma=c(2,2,10,2))
  heatmap(scRatios[c(down),], Colv=NA, margins=c(6, 9), main="Down regulated genes", col=gray.colors(4))
  op <- par(oma=c(0,0,0,0))
  legend("topleft",title="FoldChange", c("High","Medium","Low","Lowest"), col=c(gray.colors(4)),lwd=c(10))
dev.off()

# ------------------------ Create overlap PP ... IL Heatmap  ------------------------ #

PP_LN <- unlist(lapply(sharedPP_LN, function(x){which(grepl(x, rownames(ratios)))}))      # ONLY the down regulated ones !!!! because the DOWN overwrite the UP model

postscript(file="SharedPP_LN_4_tissues.eps", width=10, height=10)
  heatmap(scRatios[c(PP_LN),], margins=c(6, 9), main="Shared between LN and PP (Down regulated)", col=gray.colors(4), sub="")
  legend("topleft",title="FoldChange", c("High","Medium","Low","Lowest"), col=c(gray.colors(4)),lwd=c(10))
dev.off()

postscript(file="SharedPP_LN_4_tissues_NOTREE.eps", width=10, height=15)
  op <- par(oma=c(2,2,10,2))
  heatmap(scRatios[c(PP_LN),], Colv=NA, margins=c(6, 9), main="Shared between LN and PP (Down regulated)", col=gray.colors(4), sub="")
  op <- par(oma=c(0,0,0,0))
  legend("topleft",title="FoldChange", c("High","Medium","Low","Lowest"), col=c(gray.colors(4)),lwd=c(10))
dev.off()

# ------------------------ Create the Treatment versus Control heatmap  ------------------------ #

TvsC <- cbind(apply(counts[c(up,down),5:8],1,mean), apply(counts[c(up,down),1:4],1,mean))
colnames(TvsC) <- c("Treatment", "Control")
fullNames <- rownames(TvsC)
rownames(TvsC) <- unlist(lapply(strsplit(fullNames, ":"),"[", 2))                                   # Remove chromosome & Location
reNameToENSSSCG <- which(rownames(TvsC) == "NA")                                                    # NA gene names
rownames(TvsC)[reNameToENSSSCG] <- unlist(lapply(strsplit(fullNames[reNameToENSSSCG], ":"),"[", 1)) # Use the ENSSSCG ID in plots

TvsC <- t(apply(TvsC,1,function(x){ return(x/max(abs(x))) }))

postscript(file="Treatment_vs_Control.eps", width= 10, height = 10)
  heatmap(TvsC, scale="none", Rowv=NA, Colv=NA, margins=c(15,10), col=c(gray.colors(10)), main="Treatment vs Control")
dev.off()

Jennymatrix <- rbind(c("IGLC", 0.86, 0.61438169, 1.05031759),
c("IGLV-4", 0.39, 3.76900462, 0.25368569),
c("IGLV-12", 1.08100798, 0.90863049, 1.10661884),
c("IGKC", 0.98209508, 0.69016745, 1.24069134),
c("IL4I1", 0.83520808, 0.83520808, 0.83643003),
c("PTPRCAP", 0.54646822, 1.50614902, 1.22101453),
c("CCL17", 0.94367609, 1.06694952, 1.62117059))

colfunc <- colorRampPalette(c("green", rgb(0.9, 0.9, 0.9), "red"))

postscript(file="MyTry.jenny.eps", width=11, height=10, family = "Arial", horizontal=FALSE)

  image(t(apply(Jennymatrix[,2:4],2,as.numeric)), col=colfunc(7), breaks= c(0, 0.25, 0.5,0.9,1.1, 2, 8, 100), xaxt='n', yaxt='n')
  box()
  axis(1, at=c(0, 0.5, 1), c("Day 14","Day 35", "Day 56"))
  axis(2, at=seq(0, 1, 1/(nrow(Jennymatrix)-1)), Jennymatrix[,1],las=2)
  grid(3,7,lty=1,col="white")
  legend("topleft", title="Ratios", c("0.25 - 0.50", "0.50 - 0.90","0.90 - 1.10","1.10 - 2.00", "2.00 - 8.00"), fill=colfunc(7)[-c(1,7)])

dev.off()