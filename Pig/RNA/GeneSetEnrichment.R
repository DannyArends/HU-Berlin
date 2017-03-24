library(fgsea)
library(biomaRt)
setwd("D:/Edrive/Pig/RNA/Zink_mRNASequencing")


#mRNA <- mRNA[-which(duplicated(mRNA[,1])),]
mRNA <- read.csv("differentialExpressed.txt", sep="\t", header = TRUE)
mRNA <- mRNA[order(mRNA$pvalue), ]
write.table(mRNA[which(mRNA$pvalue < 0.05), ], "DifferentialExpression.txt", sep = "\t",row.names=FALSE)

ensembl = useMart("ensembl", "sscrofa_gene_ensembl")
conversionBioMart  <- getBM(attributes = c("ensembl_gene_id", "entrezgene"), filter="ensembl_gene_id", values=mRNA[,"ensembl_gene_id"], mart = ensembl)
conversionViaHuman <- read.csv("conv_2F28242CF0081487844190251.txt", row.names=NULL,sep="\t", colClasses="character")

mRNA <- cbind(mRNA, ENTREZ = NA)
for(x in 1:nrow(mRNA)){
  if(is.na(mRNA[x, "ENTREZ"])){
    idxConv <- which(conversionBioMart[,1] == mRNA[x,"ensembl_gene_id"])
    if(length(idxConv) == 1 ){
      mRNA[x, "ENTREZ"] <- conversionBioMart[idxConv, 2]
    }
  }
  if(is.na(mRNA[x, "ENTREZ"])){
    idxConv2 <- which(conversionViaHuman[,1] == mRNA[x,"hsapiens_homolog_ensembl_gene"])
    if(length(idxConv2) == 1){
      mRNA[x, "ENTREZ"] <- conversionViaHuman[idxConv2, 2]
    }
  }
}

filteredMRNA <- mRNA

#filteredMRNA <- filteredMRNA[-which(is.na(filteredMRNA[,"ENTREZ"])),]
#filteredMRNA <- filteredMRNA[-which(duplicated(filteredMRNA[,"ensembl_gene_id"])),]

rankings <- filteredMRNA[,"pvalue"]
names(rankings) <- filteredMRNA[,"ENTREZ"]

pathways <- reactomePathways(names(rankings))

set.seed(1)
fgseaRes <- fgsea(pathways, rankings, nperm = 10000000, maxSize=500, minSize=100)
fgseaRes <- fgseaRes[order(padj), ]

topPathwaysUp <- unlist(head(fgseaRes[NES > 0], n = 8)[,"pathway"])
#topPathwaysDown <- unlist(head(fgseaRes[NES < 0], n = 10)[,"pathway"])
topPathways <- c(unlist(topPathwaysUp[-which(duplicated(topPathwaysUp))])) #, topPathwaysDown)

#png("GSEAtable_top10up_top10down.png", width=1024, height=800)
plotGseaTable(pathways[topPathways], rankings, fgseaRes, gseaParam = 0.5)


#cat(mRNA[,"hsapiens_homolog_ensembl_gene"], sep="\n", file="ensemblID.txt")

conversion2 <- read.csv("conv_2F28242CF0081487844190251.txt", row.names=NULL,sep="\t", colClasses="character")


for(x in 1:nrow(mRNA)){
  if(is.na(mRNA[x, "ENTREZ"])){
    idxConv2 <- which(conversion2[,1] == mRNA[x,"hsapiens_homolog_ensembl_gene"])
    if(length(idxConv2) == 1){
      mRNA[x, "ENTREZ"] <- conversion2[idxConv2, 2]
    }
  }
  if(is.na(mRNA[x, "ENTREZ"])){
    idxConv <- which(conversionTable[,1] == mRNA[x,"ensembl_gene_id"])
    if(length(idxConv) == 1 ){
      mRNA[x, "ENTREZ"] <- conversionTable[idxConv, 2]
    }
  }
}

isNA <- which(is.na(mRNA[,"ENTREZ"]))

length(isNA)
mRNA <- mRNA[-isNA,]
mRNA <- mRNA[-which(duplicated(mRNA[,"ENTREZ"])),]

rankings <- mRNA[,"pvalue"]
names(rankings) <- mRNA[,"ENTREZ"]

pathways <- reactomePathways(unique(names(rankings)))

set.seed(1)
fgseaRes <- fgsea(pathways, rankings, nperm = 100000, maxSize=500, minSize=20)
fgseaRes <- fgseaRes[order(padj), ]
fgseaRes <- fgseaRes[padj < 0.05, ]
fgseaRes

topPathwaysUp <- unlist(head(fgseaRes[ES > 0], n = 10)[,"pathway"])
topPathwaysDown <- unlist(head(fgseaRes[ES < 0], n = 10)[,"pathway"])
topPathways <- c(topPathwaysUp, topPathwaysDown)

#png("GSEAtable_top10up_top10down.png", width=1024, height=800)
plotGseaTable(pathways[topPathways], rankings, fgseaRes, gseaParam = 0.5)
#dev.off()

for(x in 1:length(topPathways)){
  png(paste0(x, "_", gsub("/", "", topPathways[x]), ".png"), width=1024, height=800)
  plotEnrichment(pathways[[topPathways[x]]], rankings) + labs(title=topPathways[x])
  dev.off()
}
