# Analysis of RNA seq data for Filip Larsberg
#library("biomaRt")

setwd("D:/Edrive/Pig/RNA/S906 wholeRNAsequencing/03_exprdata")
raw_counts <- read.table("raw_counts.tsv", sep = "\t",header=TRUE, row.names=1)
normalized_counts <- read.table("normalized_counts.tsv", sep = "\t",header=TRUE, row.names=1)

notExpressed <- which(apply(normalized_counts, 1, mean) < 3)
normalized_counts <- round(normalized_counts[-notExpressed,], 1)

both <- rownames(sug.ev)[which(rownames(sug.ev) %in% rownames(sug.evuv))]
p.ev[both,]
p.evuv[both,]

#mart <- useMart(biomart = "ensembl", dataset = "sscrofa_gene_ensembl")
#annot <- getBM(attributes = c("ensembl_gene_id", 
                              #"chromosome_name", "start_position", "end_position", "strand",
                              #"external_gene_name", "wikigene_description"), filters = "ensembl_gene_id", values = rownames(normalized_counts), mart = mart)
# Write out the annotation to file, so we dont have to do the download a second time
#write.table(annot, "annotation.txt", sep="\t", quote=FALSE, row.names=FALSE)

annot <- read.table("annotation.txt", sep="\t", header=TRUE, quote = "", colClasses = "character")

ctrl <- c(1, 4, 7, 10)
ef <- c(2, 5, 8, 11)
efuv <- c(3, 6, 9, 12)

p.ev <- t(apply(normalized_counts,1,function(x){
  return(c(mean(x[ctrl]), round(sd(x[ctrl]),1), mean(x[ef]), round(sd(x[ef]),1), mean(x[ef])/mean(x[ctrl]), log2(mean(x[ef])/mean(x[ctrl])), wilcox.test(x[ctrl], x[ef])$p.value))
}))
colnames(p.ev) <- c("mean(ctrl)", "sd(ctrl)", "mean(ev)", "sd(ev)", "ratio", "log2(ratio)", "p.value")

p.evuv <- t(apply(normalized_counts,1,function(x){
  return(c(mean(x[ctrl]), round(sd(x[ctrl]),1), mean(x[efuv]), round(sd(x[efuv]),1), mean(x[efuv])/mean(x[ctrl]), log2(mean(x[efuv])/mean(x[ctrl])), wilcox.test(x[ctrl], x[efuv])$p.value))
}))
colnames(p.evuv) <- c("mean(ctrl)", "sd(ctrl)", "mean(evuv)", "sd(evuv)", "ratio", "log2(ratio)", "p.value")

diff.ev <- names(which(p.ev[,"p.value"] < 0.05))
diff.evuv <- names(which(p.evuv[,"p.value"] < 0.05))

sug.ev <- p.ev[diff.ev,]
sug.ev <- sug.ev[sort(sug.ev[,"log2(ratio)"], index.return=TRUE)$ix,]
ensgid <- rownames(sug.ev)
annot.ev <- lapply(ensgid, function(x){
  return(annot[which(annot[,1] == x)[1],])
})
sug.ev <- cbind(chr = unlist(lapply(annot.ev, "[", 2)), sug.ev)
sug.ev <- cbind(start = unlist(lapply(annot.ev, "[", 3)), sug.ev)
sug.ev <- cbind(stop = unlist(lapply(annot.ev, "[", 4)), sug.ev)
sug.ev <- cbind(strand = unlist(lapply(annot.ev, "[", 5)), sug.ev)
sug.ev <- cbind(description = unlist(lapply(annot.ev, "[", 7)), sug.ev)
sug.ev <- cbind(genename = unlist(lapply(annot.ev, "[", 6)), sug.ev)
rownames(sug.ev) <- ensgid
write.table(sug.ev, "ev_vs_ctrl_0.05.txt", sep="\t", quote=FALSE)

sug.evuv <- p.evuv[diff.evuv,]
sug.evuv <- sug.evuv[sort(sug.evuv[,"log2(ratio)"], index.return=TRUE)$ix,]
ensgid <- rownames(sug.evuv)
annot.evuv <- lapply(ensgid, function(x){
  return(annot[which(annot[,1] == x)[1],])
})

sug.evuv <- cbind(chr = unlist(lapply(annot.evuv, "[", 2)), sug.evuv)
sug.evuv <- cbind(start = unlist(lapply(annot.evuv, "[", 3)), sug.evuv)
sug.evuv <- cbind(stop = unlist(lapply(annot.evuv, "[", 4)), sug.evuv)
sug.evuv <- cbind(strand = unlist(lapply(annot.evuv, "[", 5)), sug.evuv)
sug.evuv <- cbind(description = unlist(lapply(annot.evuv, "[", 7)), sug.evuv)
sug.evuv <- cbind(genename = unlist(lapply(annot.evuv, "[", 6)), sug.evuv)
rownames(sug.evuv) <- ensgid
write.table(sug.evuv, "evuv_vs_ctrl_0.05.txt", sep="\t", quote=FALSE)

mir_ev <- cbind(annot[grep("ssc-mir", annot[,"external_gene_name"]),], p.ev[annot[grep("ssc-mir", annot[,"external_gene_name"]),1],])
write.table(mir_ev, "miRNA_ev.txt", sep="\t", quote=FALSE)

mir_evuv <- cbind(annot[grep("ssc-mir", annot[,"external_gene_name"]),], p.evuv[annot[grep("ssc-mir", annot[,"external_gene_name"]),1],])
write.table(mir_evuv, "miRNA_evuv.txt", sep="\t", quote=FALSE)
