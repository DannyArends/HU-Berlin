# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

library(biomaRt)

setwd("E:/Mouse/RNA/FV3/Analysis_Helge")
processed <- read.csv("FV3_ttest_final.txt", skip=56, sep="\t", header = TRUE, na.strings="\\N")
expressed <- processed[which(processed[,"p_HFD_liver"] < 2),]
expressed <- expressed[which(expressed[,"blast_Ensembl_Gene_Id"] != ""),]
genes <- unique(as.character(na.omit(expressed[,"blast_Ensembl_Gene_Id"])))

mart      <- useMart("ensembl", "mmusculus_gene_ensembl")
allgenes  <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),  filter="ensembl_gene_id", values=genes, mart = mart)

myexp <- NULL
for(x in 1:nrow(expressed)){
  ensembl_gene_id <- as.character(expressed[x,"blast_Ensembl_Gene_Id"])
  mgi_symbol <- allgenes[which(allgenes[,"ensembl_gene_id"] == ensembl_gene_id),"mgi_symbol"]
  pvalue <- expressed[x,"p_HFD_liver"]
  foldchange <- expressed[x,"fch_HFD_liver"]
  if(length(mgi_symbol) == 0) mgi_symbol <- ""
  xrow <- c(ensembl_gene_id, mgi_symbol, pvalue, foldchange)
  cat(x, length(xrow),"\n")
  myexp <- rbind(myexp, xrow)
}

colnames(myexp) <- c("ensembl_gene_id", "mgi_symbol", "pvalue", "foldchange")

myexp <- myexp[sort(as.numeric(myexp[,"pvalue"]), index.return=TRUE)$ix,]

write.table(myexp[!duplicated(myexp[,"ensembl_gene_id"]), ], file="Helge_ensembl_expressed_3e-5.txt", sep ="\t", quote=FALSE, row.names=FALSE)
