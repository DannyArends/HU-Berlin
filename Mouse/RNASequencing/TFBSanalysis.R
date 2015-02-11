# RNA Seq - Expression data analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015

library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
mdata <- read.table("Analysis/TFBSmatrix.txt", sep = "\t", header=TRUE)
genes <- unique(mdata[mdata[,"decision.majority.IGV"] == "B6N","ensembl_gene_id"])

chromosomal.loc <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene") [genes[-1]]
promoter.gene   <- getPromoterSeq(chromosomal.loc, Mmusculus, upstream=2000, downstream=1000)

TFBSjasparALL <- query(MotifDb, "Mmusculus")

output <- matrix(NA,  length(names(TFBSjasparALL)),length(genes), dimnames=list(names(TFBSjasparALL), genes))

for(tfbs in 1:length(TFBSjasparALL)){
  seqLogo(TFBSjasparALL[[tfbs]])
  TFBSjaspar <- round(100 * TFBSjasparALL[[tfbs]])
  for(gene in names(promoter.gene)){
    hits <- matchPWM(TFBSjaspar, unlist(promoter.gene[gene])[[1]], "90%")
    output[names(TFBSjasparALL)[tfbs], gene] <- length(hits)
  }
}
write.table(output, file="Analysis/TFBSmatches_B6N.txt",sep="\t")


TBFsitesB6N <-  read.table("Analysis/TFBSmatches_B6N.txt", sep = "\t", header=TRUE)
TBFsitesBFMI <-  read.table("Analysis/TFBSmatches_BFMI.txt", sep = "\t", header=TRUE)

overview <- cbind(B6N = apply(TBFsitesB6N,1,function(x){sum(x != 0,na.rm=TRUE)/length(na.omit(x))}), BFMI = apply(TBFsitesBFMI,1,function(x){sum(x != 0,na.rm=TRUE)/length(na.omit(x))}))
write.table(round(overview, d = 3), file="Analysis/TFBSoverview.txt",sep="\t")


