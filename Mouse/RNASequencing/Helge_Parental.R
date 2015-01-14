# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

library(biomaRt)

setwd("E:/Mouse/RNA/FV3/Analysis_Helge")
processed <- read.csv("FV3_ttest_final.txt", skip=56, sep="\t", header = TRUE, na.strings="\\N")
expressed <- processed[which(processed[,"p_HFD_fat"] < 0.0000003),]

cat(unique(as.character(na.omit(expressed[,"blast_Ensembl_Gene_Id"]))),sep="\n", file="Helge_ensembl_expressed_3e-7.txt")

#mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
#res.biomart <- getBM(attributes = c("illumina_mousewg_6_v1", "mgi_id","mgi_symbol", "chromosome_name", "start_position", "end_position", "strand"), filters = "illumina_mousewg_6_v1", values = expressed[,"ID"], mart = mart)
