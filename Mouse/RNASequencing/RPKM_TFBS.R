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

#RPKM <- read.table(file="Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep="\t")

#SUBSET <- RPKM[which(RPKM[,"Ratio_F1_Scale"] > 1.2 | RPKM[,"Ratio_F1_Scale"] < -1.2),]
#SUBSET <- SUBSET[which(SUBSET[,"tTest_F1"] < 0.1),]
#SUBSET <- SUBSET[which(SUBSET[,"Mean.BFMI860.12xB6N.L"] > 0.5 | SUBSET[,"Mean.B6NxBFMI860.12.L"] > 0.5),]

#genes <- SUBSET[,"ensembl_gene_id"]
genes <- as.character(unlist(read.table("S2_ENS.txt")))
#genes <- genes[-which(!genes %in% names(transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene")))]

chromosomal.loc <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene") [genes]
promoter.gene   <- getPromoterSeq(chromosomal.loc, Mmusculus, upstream=2000, downstream=1000)

TFBSjasparALL <- query(MotifDb, "Mmusculus")

output <- matrix(NA,  length(names(TFBSjasparALL)),length(genes), dimnames=list(names(TFBSjasparALL), genes))

estrogen <- c("Mmusculus-jolma2013-Esrra-2", "Mmusculus-UniPROBE-Esrra.UP00079", "Mmusculus-JASPAR_CORE-Esrrb-MA0141.1", "Mmusculus-JASPAR_2014-Esrrb-MA0141.2")

for(tfbs in 1:length(TFBSjasparALL)){
  seqLogo(TFBSjasparALL[[tfbs]])
  TFBSjaspar <- round(100 * TFBSjasparALL[[tfbs]])
  for(gene in names(promoter.gene)){
    hits <- matchPWM(TFBSjaspar, unlist(promoter.gene[gene])[[1]], "90%")
    output[names(TFBSjasparALL)[tfbs], gene] <- length(hits)
  }
}


library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes

fullnames <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), filters = c("ensembl_gene_id"), values = colnames(output), mart = bio.mart)

output2 <- output[,match(fullnames[,"ensembl_gene_id"], colnames(output))]

write.table(output, file="Analysis/TFBSmatches_RPKM.txt",sep="\t")
write.table(rbind(fullnames[,"mgi_symbol"], output2), file="Analysis/TFBSmatches_S2_ENS_ALL.txt",sep="\t",quote=FALSE)

weights.ADR3 <- rbind(c(1,0,1,1,0,1,0.25,0.25,0.25,1,0,1,1,0,1),        #A
                      c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,0),        #C
                      c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0),        #G
                      c(0,0,0,0,0,0,0.25,0.25,0.25,0,0,0,0,0,0))        #T
rownames(weights.ADR3) <- c("A","C","G","T")
weights.ARE2 <- rbind(c(1,0,0,1,0,0,0.25,0.25,0.25,1,0,1,1,0,0),        #A
                      c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,1),        #C
                      c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0),        #G
                      c(0,0,1,0,0,1,0.25,0.25,0.25,0,0,0,0,0,0))        #T
rownames(weights.ARE2) <- c("A","C","G","T")
weights.ARE  <- rbind(c(1,0,1,1,0,1,0.25,0.25,0.25,0,0,0,1,0,0),        #A
                      c(0,0,0,0,1,0,0.25,0.25,0.25,1,0,0,0,1,0),        #C
                      c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0),        #G
                      c(0,0,0,0,0,0,0.25,0.25,0.25,0,0,1,0,0,1))        #T
rownames(weights.ARE) <- c("A","C","G","T")
weights.IR3  <- rbind(c(1,0,1,1,0,1,0.25,0.25,0.25,0,0,0,0,0,0),        #A
                      c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,0),        #C
                      c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,1),        #G
                      c(0,0,0,0,0,0,0.25,0.25,0.25,1,0,1,1,0,0))        #T
rownames(weights.IR3) <- c("A","C","G","T")

wmatrices <- list(weights.ADR3, weights.ARE2, weights.ARE, weights.IR3)
names(wmatrices) <- c("ADR3","ARE2","ARE","IR3")

output <- NULL
for(x in 1:length(wmatrices)){
  seqLogo(wmatrices[[x]])
  TFBSjaspar <- round(100 *wmatrices[[x]])
  allhits <- NULL
  for(gene in names(promoter.gene)){
    hits <- matchPWM(TFBSjaspar, unlist(promoter.gene[gene])[[1]], "90%")
    allhits <- c(allhits, length(hits))
  }
  output <- cbind(output, allhits)
}
rownames(output) <- names(promoter.gene)
colnames(output) <- names(wmatrices)

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes

fullnames <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), filters = c("ensembl_gene_id"), values = rownames(output), mart = bio.mart)

output <- output[match(fullnames[,"ensembl_gene_id"], colnames(output)),]

write.table(cbind(fullnames[,"mgi_symbol"], output), file="Analysis/TFBSmatches_Estrogen.txt",sep="\t",quote=FALSE)


fullnames <- getBM(attributes = c("ensembl_gene_id", "ggallus_homolog_ensembl_gene"), filters = c("ensembl_gene_id"), values = "ENSMUSG00000059201", mart = bio.mart)


