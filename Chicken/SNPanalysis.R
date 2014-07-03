# Analysis of Chicken 600K SNP chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

setwd("E:/Chicken/DNA/600KSNPChip/")

calldata <- read.table("GenotypeData.txt", na.strings="-1")
chipids <- read.table("Annotation/originalChipIDs.txt", sep="\t", header=TRUE)

dendrogram <- as.dendrogram(hclust(dist(t(calldata))))                       # Dendrogram of all ALL lines
plot(dendrogram, main="Dendrogram All Chicken", xlab="", center=TRUE)

indPop1 <- paste0("X", chipids[which(chipids[,"popPk"]=="IBNHxx"),"ID_Chip"])
indPop2 <- paste0("X", chipids[which(chipids[,"popPk"]=="IB77xx"),"ID_Chip"])

nGenoPop1 <- apply(calldata[,indPop1],1,function(x){return(length(na.omit(unique(x))))})
nGenoPop2 <- apply(calldata[,indPop2],1,function(x){return(length(na.omit(unique(x))))})

mightBeDifferent <- which((nGenoPop1 + nGenoPop2) > 1 & (nGenoPop1 + nGenoPop2) < 4)

calldata <- calldata[mightBeDifferent,]
genoInP1 <- apply(calldata[,indPop1], 1, function(x){return(as.numeric(na.omit(unique(x))))})
genoInP2 <- apply(calldata[,indPop2], 1, function(x){return(as.numeric(na.omit(unique(x))))})

different <- rep(FALSE, nrow(calldata))
for(x in 1:nrow(calldata)){
  if(length(genoInP1[[x]])==1) different[x] = !(genoInP1[[x]] %in% genoInP2[[x]])
  if(length(genoInP2[[x]])==1) different[x] = !(genoInP2[[x]] %in% genoInP1[[x]])
  if(x %% 10000 == 0)cat("Done",x,"/",length(mightBeDifferent),"\n")
}

calldata <- calldata[which(different),]
write.table(calldata[,c(indPop1,indPop2)], file="Analysis/DifferentialGenotypes.txt", sep="\t")

arrayAnnotation <- read.table("Annotation/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", header=TRUE, colClasses="character", na.string="---")

annotSubset <- which(arrayAnnotation[,"Probe.Set.ID"] %in% rownames(calldata))
arrayAnnotation <- arrayAnnotation[annotSubset,]

affyRelations <- c("splice-site", "intron", "nonsense", "missense", "synon", "frameshift", "aa-indel", "CDS", "UTR-3", "UTR-5", "exon", "upstream", "downstream")

strcount <- function(x, pattern, split){ unlist(lapply(strsplit(x, split), function(z) na.omit(length(grep(pattern, z))))) }

relationCount <- NULL
for(relation in affyRelations){
  inSNP <- grep(relation, arrayAnnotation[,"Associated.Gene"])
  cnt <- 0
  for(line in arrayAnnotation[inSNP,"Associated.Gene"]){
    cnt = cnt + strcount(line, relation, " // ")
  }
  relationCount <- rbind(relationCount, c(length(unique(inSNP)), cnt))
}
rownames(relationCount) <- affyRelations
colnames(relationCount) <- c("SNPs","Genes")
