#
# Post analysis of ZInk data
#
setwd("E:/Pig/RNA/Zink_mRNASequencing")
rawreads <- read.table("RawReads.txt", sep="\t", check.names=FALSE)
rpkm <- read.table("RPKMnorm.txt", sep="\t", check.names=FALSE)
colnames(rpkm) <- unlist(lapply(strsplit(colnames(rpkm),".",fixed=TRUE), "[",1))

# Groups
medium <- c("91", "112", "115"); high <- c("92", "113", "116")

rpkm[rpkm == 0] <- NA
rpkm <- log2(rpkm)
rpkm[is.na(rpkm)] <- 0

cat("Before:", nrow(rpkm), "genes\n")
rpkm <- rpkm[-which(apply(rpkm[,medium],1,function(x){all(x < 1)})),]
rpkm <- rpkm[-which(apply(rpkm[,high],1,function(x){all(x < 1)})),]

cat("After:", nrow(rpkm), "genes\n")

# How do the samples cluster
clusters <- hclust(dist(t(rpkm)))
plot(clusters)

# Histogram of variation across the genes
variance <- apply(rpkm,1, var)
variance[variance == 0] <- NA
logvar <- abs(log(variance))
hist(logvar)

# Most variable genes across the population
mostvariable <- rownames(rpkm)[which(logvar > 7)]


library(biomaRt)
mart = useMart(biomart="ensembl", dataset="sscrofa_gene_ensembl")
genenames  <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "hsapiens_homolog_ensembl_gene"), 
                    filter = "ensembl_gene_id", values = mostvariable, mart = mart)

write.table(genenames, sep="\t", file="mostVariableGenes.txt", quote=FALSE, row.names=FALSE)                    
cat(unique(genenames[,3]), sep="\n", file="humanOrthologs.txt")


pvals <- t(apply(cbind(rpkm[,medium],rpkm[,high]),1,function(x){
  return(t(c(round(mean(x[1:3]),2), round(mean(x[4:6]),2), t.test(x[1:3], x[4:6])$p.value)))
}))

genenames  <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description", "hsapiens_homolog_ensembl_gene"), 
                    filter = "ensembl_gene_id", values = rownames(pvals), mart = mart)

annotatedDE <- cbind(genenames[,1], pvals[genenames[,1],], genenames[,2:4])
colnames(annotatedDE) <- c("ensembl_gene_id", "mean(control)", "mean(treatment)", "pvalue", "external_gene_name", "description", "hsapiens_homolog_ensembl_gene")

# Sort at pvalue
annotatedDE <- annotatedDE[sort(annotatedDE[,"pvalue"], index.return=TRUE)$ix,]
p001 <- annotatedDE[which(annotatedDE[,"pvalue"] < 0.01), ]
write.table(annotatedDE, sep="\t", file="differentialExpressed.txt", quote=FALSE, row.names=FALSE)
write.table(p001, sep="\t", file="differentialExpressed0p01.txt", quote=FALSE, row.names=FALSE)

cat(unique(p001[,"hsapiens_homolog_ensembl_gene"]), sep="\n", file="humanOrthologsDE.txt")

library(biomaRt)
mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

humanBM  <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"), 
                    filter = "ensembl_gene_id", values = unique(p001[,"hsapiens_homolog_ensembl_gene"]), mart = mart)

p001 <- cbind(p001, "human_gene_name" = NA, "human_gene_description" = NA)
for(x in 1:nrow(p001)){
  idx <- which(humanBM[,1] == p001[x,"hsapiens_homolog_ensembl_gene"])
  if(length(idx) == 1){
    p001[x,"human_gene_name"] <- humanBM[idx,2]
    p001[x,"human_gene_description"] <- humanBM[idx,3]
  }
}
write.table(p001, sep="\t", file="differentialExpressed0p01humanAnnot.txt", quote=FALSE, row.names=FALSE)
