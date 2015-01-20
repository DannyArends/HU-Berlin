# Analysis of FV3 RNA expression data collected using Illumina microarrays
# Finding differences between BFMI and B6 in Liver under FF
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

setwd("E:/Mouse/RNA/FV3")

alldata  <- read.table("Analysis/geneexpression.txt", sep="\t", header=TRUE, check.names=FALSE)
arrays   <- read.table("Annotation/arrays.txt", sep="\t", header=TRUE, row.names=1)                                   # Arrays annotation

FFliver <- arrays[which(arrays[,"tissue"] == "liver" & arrays[,"diet"] == "FF"), ]
B6 <- which(FFliver[,"line"]=="B6")
BFMI <- which(FFliver[,"line"]=="BFMI860")

pvalues <- t(apply(alldata[,rownames(FFliver)], 1, function(x){
  mB6 <- mean(x[B6])
  mBFMI <- mean(x[BFMI])
  if(mB6 > mBFMI){
    scR <- -(mB6/mBFMI)
  }else{
    scR <- mBFMI / mB6
  }
  c(round(mB6,digits=2), round(mBFMI, digits=2), t.test(x[B6],x[BFMI])$p.value, mBFMI / mB6, scR)
}))
colnames(pvalues) <- c("B6", "BFMI", "Pvalue", "Ratio", "Ratio_Scale")
pvalues <- cbind(ProbeName = as.character(alldata[,"ProbeName"]), ensembl_gene_id = as.character(alldata[,"ensembl_gene_id"]), mgi_symbol = as.character(alldata[,"mgi_symbol"]), pvalues)

top400 <- pvalues[sort(as.numeric(pvalues[,"Pvalue"]),index.return=TRUE)$ix[1:400],]
top400 <- top400[!duplicated(top400[,"ensembl_gene_id"]), ]
write.table(top400, file="Analysis/Helge_ensembl_top400.txt", sep="\t", row.names = FALSE, quote=FALSE)


pvalues[which(pvalues[,"ProbeName"] == "1090142"),]

alldata[which(alldata[,"ProbeName"] == "100070026"),]
