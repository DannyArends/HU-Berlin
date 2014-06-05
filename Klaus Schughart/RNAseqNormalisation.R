# RNAseqNormalisation.R
#
# copyright (c) 2014-2020 - Junxi Wang, Klaus Schughart and Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Script for generating the GOOD PCAs (171213) with DESeq

library(DESeq2)         # Load the DESeq2 package
setwd("d:/Klaus_Mouse") # Dort sind die Daten

# Read data
all_mouse_counts  <- read.table("all_mouse_counts.txt", header=TRUE)

# Subset the data, to only contain the values
rawdata           <- all_mouse_counts[, 8:ncol(all_mouse_counts)] # Subset data starting from colum 8

# Weil einige Gen-Namen gleich sind, nehme ich die Nummern der Reihen als Marker für jede Reihe.
rownames(rawdata) <- all_mouse_counts[,1]

names(rawdata)  # Check the column names
dim(rawdata)    # Check the dimensions

# Damit gabe es kein Null, deshalb wird kein Log2FolgChange-Wert als NA.
# Add 1 to the expression values so we can do log(exp), without erroring on 0
neudat <- rawdata + 1 

names(neudat)   # Check the column names
dim(neudat)     # Check the dimensions

summary(neudat)
boxplot(neudat)

# Define the clusters
cluster <- c("Count.B6_Mock_d1", "Count.B6_Mock_d3", "Count.B6_Inf_d1", "Count.B6_Inf_d3", "Count.B6_Inf_d5", "Count.B6_Inf_d8",
             "Count.B6_Inf_d14", "Count.D2_Mock_d1","Count.D2_Mock_d3","Count.D2_Inf_d1","Count.D2_Inf_d3","Count.D2_Inf_d5")

conditions <- NULL                # Danny: better to use NULL not c()
for(i in 1:length(cluster)){
  conditions <- c(conditions, cluster[i], cluster[i], cluster[i])
}

# Create the DEseq2 Dataset structure
columnData <- data.frame(condition=factor(conditions))
dsDat      <- DESeqDataSetFromMatrix(countData = neudat, colData = columnData, design = ~condition)

# Transformation and PCA-Analyse
rld <- rlogTransformation(dsDat, blind=TRUE)
print(plotPCA(rld))

# Get the logTransformed data out
datafromRld <- as.data.frame(assays(rld)[[1]])                                   # Get the data from the RLD object, as a dataframe
datafromRld <- cbind(as.character(all_mouse_counts[,"Gene.ID"]), datafromRld)    # Add an additional column for the gene.ID
rownames(datafromRld) <- rownames(rawdata)                                       # Set the rownames to the rownames of rawdata/neudata
colnames(datafromRld) <- c("Gene.ID", names(rawdata))                            # Add column names, such that we can id the samples

# Write it twice ?
write.table(datafromRld, file="RNAseq_norm1.txt")                                # store normalized data as plain text
write.csv2(as.data.frame(datafromRld), "RNAseq_norm1.csv")                       # Store the normalized data as a csv
