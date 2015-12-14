setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross/Analysis")
expdata <- read.table(file=paste0("RPKM_norm_log_stats.txt"), sep="\t", header=TRUE)
rownames(expdata) <- expdata[,"ensembl_gene_id"]

## Quadriceps
rr <- which(expdata[,"pValues_Q"] < 0.05 & (expdata[,"ratios_Q"] < 0.83 | expdata[,"ratios_Q"] > 1.2))
cc <- grep("_Q", colnames(expdata))
#cc <- cc[cc < 35][1:6]
cc <- cc[cc < 35]

correlations <- NULL
for(x in rr){
  correlations <- rbind(correlations, cor(as.numeric(expdata[x,cc]), t(expdata[rr,cc])))
}
rownames(correlations) <- colnames(correlations)
write.table(correlations, "Correlations_Q.txt", sep="\t")
heatmap(correlations, main="Quadriceps")

## Liver
rr <- which(expdata[,"pValues_L"] < 0.05 & (expdata[,"ratios_L"] < 0.83 | expdata[,"ratios_L"] > 1.2))
cc <- grep("_L", colnames(expdata))
cc <- cc[cc < 35]

correlations <- NULL
for(x in rr){
  correlations <- rbind(correlations, cor(as.numeric(expdata[x,cc]), t(expdata[rr,cc])))
}
rownames(correlations) <- colnames(correlations)
write.table(correlations, "Correlations_L.txt", sep="\t")
heatmap(correlations, main="Liver")

## Gonadal Fat
rr <- which(expdata[,"pValues_G"] < 0.05 & (expdata[,"ratios_G"] < 0.83 | expdata[,"ratios_G"] > 1.2))
cc <- grep("_G", colnames(expdata))
cc <- cc[cc < 35]

correlations <- NULL
for(x in rr){
  correlations <- rbind(correlations, cor(as.numeric(expdata[x,cc]), t(expdata[rr,cc])))
}
rownames(correlations) <- colnames(correlations)
write.table(correlations, "Correlations_G.txt", sep="\t")
heatmap(correlations, main="Gonadal Fat")

cols <- as.numeric(grepl("^B6N_", colnames(expdata)[cc]))
cols <- cols + (as.numeric(grepl("^BFMI860.12_", colnames(expdata)[cc])) * 2)
cols <- cols + (as.numeric(grepl("^BFMI860.12xB6N", colnames(expdata)[cc])) * 3)
cols <- cols + (as.numeric(grepl("^B6NxBFMI860.12", colnames(expdata)[cc])) * 4)

plot(t(expdata[c("ENSMUSG00000065629", "ENSMUSG00000104827"), cc]), col=c("gray", "orange", "blue", "green")[cols], pch=19)
plot(t(expdata[c("ENSMUSG00000065629", "ENSMUSG00000026072"), cc]), col=c("gray", "orange", "blue", "green")[cols], pch=19)
plot(t(expdata[c("ENSMUSG00000065629", "ENSMUSG00000025981"), cc]), col=c("gray", "orange", "blue", "green")[cols], pch=19)

matBFMI <- which(grepl("^BFMI860.12xB6N", colnames(expdata)))
matB6N <- which(grepl("^B6NxBFMI860.12", colnames(expdata)))

## All tissues combined

good <- which(apply(expdata[,c(matBFMI,matB6N)],1, function(x){all(x > 0.5) })) # & expdata[,"pValues_DOC"] < 0.001)

corA <- NULL
for(x in good){
  corA <- rbind(corA, cor(as.numeric(expdata[x,c(matBFMI)]), t(expdata[good,c(matBFMI)])))
}
rownames(corA) <- colnames(corA)
rownames(corA) <- expdata[rownames(corA),"gene_name"]
colnames(corA) <- expdata[colnames(corA),"gene_name"]
write.table(corA, "matBFMI_correlation.txt",sep = "\t")

png("c_A.png", width=1024, height=1024)
op <- par(mai = c(5,5,5,5))
heatmap(corA, col=c("red", "white", "blue"), breaks=c(-1000, -0.5,0.5, 1000), scale="none", Colv=NA, Rowv=NA, main="matBFMI",cexRow = 2,cexCol = 2,margins=c(12,8))
dev.off()

corB <- NULL
for(x in good){
  corB <- rbind(corB, cor(as.numeric(expdata[x,c(matB6N)]), t(expdata[good,c(matB6N)])))
}
rownames(corB) <- colnames(corB)
rownames(corB) <- expdata[rownames(corB),"gene_name"]
colnames(corB) <- expdata[colnames(corB),"gene_name"]
write.table(corB, "matB6N_correlation.txt",sep = "\t")

png("c_B.png", width=1024, height=1024)
heatmap(corB, col=c("red", "white", "blue"), breaks=c(-1000, -0.5,0.5, 1000), scale="none", Colv=NA, Rowv=NA, main="matB6N",cexRow = 2,cexCol = 2,margins=c(12,8))
dev.off()

corDiff <- (corA - corB)
write.table(corDiff, "DIFF_correlation.txt",sep = "\t")
png("c_D.png", width=1024, height=1024)
heatmap(corDiff, main="All tissues", col=c("red", "white", "blue"), breaks=c(-1000, -0.3,0.3, 1000), scale="none", Colv=NA, Rowv=NA,cexRow = 2,cexCol = 2,margins=c(12,8))
dev.off()
