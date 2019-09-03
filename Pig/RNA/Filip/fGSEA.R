library(reactome.db)
library(org.Ss.eg.db)
library(data.table)
library(fgsea)

pig.universe <- keys(org.Ss.eg.db, "ENTREZID")

# Selecting reactome gene sets
pathways <- na.omit(select(reactome.db, keys=pig.universe, c("PATHID"),
                           keytype = 'ENTREZID'))
pathways <- split(pathways$ENTREZID, pathways$PATHID)

pathway2name <- as.data.table(na.omit(select(reactome.db, names(pathways),
                                             c("PATHNAME"), 'PATHID')))
# Remove organism prefix
pathway2name[, PATHNAME := sub("^[^:]*: ", "", PATHNAME)]
pathway2name <- setNames(pathway2name$PATHNAME, pathway2name$PATHID)
pathway2name <- iconv(pathway2name, "latin1", "ASCII", sub="")

pathway.lines <- sapply(names(pathways), function(p) {
    link <- p
    name <- paste0(p, "_", pathway2name[p])
    name <- gsub("[ ()/]+", "_", name)
    sprintf("%s\t%s\t%s", name, link, paste0(pathways[[p]], collapse="\t"))
})

SsPathway <- strsplit(pathway.lines, "\t")
names(SsPathway) <- unlist(lapply(SsPathway,"[", 1))
SsPathway <- lapply(SsPathway, function(x) x[-c(1:2)] )


setwd("D:/Edrive/Pig/RNA/S906 wholeRNAsequencing/03_exprdata/DiffExpression_based_on_normalized_counts")

ev_ctrl <- read.csv("ef_vs_ctrl.all.txt", sep = "\t")
ev_evuv <- read.csv("ef_vs_efuv.all.txt", sep = "\t")
evuv_vs_ctrl <- read.csv("efuv_vs_ctrl.all.txt", sep = "\t")

# EF versus Control

ev_ctrl <- ev_ctrl[sort(ev_ctrl[,"log2.ratio."], index.return=TRUE)$ix, ]

library(biomaRt)
bio.mart = useMart("ensembl", "sscrofa_gene_ensembl")
res.biomart <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = c("ensembl_gene_id"), values = rownames(ev_ctrl), mart = bio.mart)
new_ev_ctrl <- NULL
for(x in 1:nrow(ev_ctrl)){
  ii <- which(res.biomart[,1] == rownames(ev_ctrl)[x])
  if(length(ii) == 1){
    if(!is.na(res.biomart[ii,2])){
      new_ev_ctrl <- rbind(new_ev_ctrl, c(rownames(ev_ctrl)[x], res.biomart[ii,2], ev_ctrl[x,]))
    }
  }
}
ev_ctrl_ranks <- unlist(new_ev_ctrl[,"log2.ratio."])
names(ev_ctrl_ranks) <- new_ev_ctrl[,2]
fgseaRes <- fgsea(SsPathway, ev_ctrl_ranks, nperm=10000, maxSize=100)
fgseaResOrdered <- fgseaRes[head(order(pval), n=25)]
write.table(as.matrix(data.frame(fgseaResOrdered[which(fgseaResOrdered[,"size"] > 1),])), file="topGSEA_ef_vs_ctrl.txt",sep="\t",quote=FALSE, row.names=FALSE)


# EF versus EFuv
ev_evuv <- ev_evuv[sort(ev_evuv[,"log2.ratio."], index.return=TRUE)$ix, ]
res.biomart <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = c("ensembl_gene_id"), values = rownames(ev_evuv), mart = bio.mart)
new_ev_evuv <- NULL
for(x in 1:nrow(ev_evuv)){
  ii <- which(res.biomart[,1] == rownames(ev_evuv)[x])
  if(length(ii) == 1){
    if(!is.na(res.biomart[ii,2])){
      new_ev_evuv <- rbind(new_ev_evuv, c(rownames(ev_evuv)[x], res.biomart[ii,2], ev_evuv[x,]))
    }
  }
}
ev_evuv_ranks <- unlist(new_ev_evuv[,"log2.ratio."])
names(ev_evuv_ranks) <- new_ev_evuv[,2]
fgseaRes <- fgsea(SsPathway, ev_evuv_ranks, nperm=10000, maxSize=100)
fgseaResOrdered <- fgseaRes[head(order(pval), n=25)]
write.table(as.matrix(data.frame(fgseaResOrdered[which(fgseaResOrdered[,"size"] > 1),])), file="topGSEA_ef_vs_efuv.txt",sep="\t",quote=FALSE, row.names=FALSE)

# EFuv versus Control
evuv_vs_ctrl <- evuv_vs_ctrl[sort(evuv_vs_ctrl[,"log2.ratio."], index.return=TRUE)$ix, ]
res.biomart <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), filters = c("ensembl_gene_id"), values = rownames(evuv_vs_ctrl), mart = bio.mart)
new_evuv_vs_ctrl <- NULL
for(x in 1:nrow(evuv_vs_ctrl)){
  ii <- which(res.biomart[,1] == rownames(evuv_vs_ctrl)[x])
  if(length(ii) == 1){
    if(!is.na(res.biomart[ii,2])){
      new_evuv_vs_ctrl <- rbind(new_evuv_vs_ctrl, c(rownames(evuv_vs_ctrl)[x], res.biomart[ii,2], evuv_vs_ctrl[x,]))
    }
  }
}
evuv_vs_ctrl_ranks <- unlist(new_evuv_vs_ctrl[,"log2.ratio."])
names(evuv_vs_ctrl_ranks) <- new_evuv_vs_ctrl[,2]
fgseaRes <- fgsea(SsPathway, evuv_vs_ctrl_ranks, nperm=10000, maxSize=100)
fgseaResOrdered <- fgseaRes[head(order(pval), n=25)]
write.table(as.matrix(data.frame(fgseaResOrdered[which(fgseaResOrdered[,"size"] > 1),])), file="topGSEA_efuv_vs_ctrl.txt",sep="\t",quote=FALSE, row.names=FALSE)
