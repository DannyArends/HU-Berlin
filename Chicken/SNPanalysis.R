# Analysis of Chicken 600K SNP chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

setwd("E:/Chicken/DNA/600KSNPChip/")

calldata <- read.table("GenotypeData.txt", na.strings="-1")
chipids <- read.table("Annotation/originalChipIDs.txt", sep="\t", header=TRUE)

dendrogram <- as.dendrogram(hclust(dist(t(calldata))))                        # Dendrogram of all ALL lines
phylotree <- as.phylo(hclust(dist(t(calldata))))                              # Phylogram of all ALL lines
plot(dendrogram, type="triangle", main="Dendrogram All Chicken", xlab="", center=TRUE)
plot(phylotree, type="fan", main="Fan Phylo All Chicken")

indPop1 <- paste0("X", chipids[which(chipids[,"popPk"]=="IBNHxx"),"ID_Chip"])                       # Individuals belonging to the first strain
indPop2 <- paste0("X", chipids[which(chipids[,"popPk"]=="IB77xx"),"ID_Chip"])                       # Individuals belonging to the first strain

nGenoPop1 <- apply(calldata[,indPop1],1,function(x){return(length(na.omit(unique(x))))})            # Number of genotypes in population 1
nGenoPop2 <- apply(calldata[,indPop2],1,function(x){return(length(na.omit(unique(x))))})            # Number of genotypes in population 2

mightBeDifferent <- which((nGenoPop1 + nGenoPop2) > 1 & (nGenoPop1 + nGenoPop2) < 4)                # Might be different between the two RI strains

calldata <- calldata[mightBeDifferent,]
genoInP1 <- apply(calldata[,indPop1], 1, function(x){return(as.numeric(na.omit(unique(x))))})       # Genotypes in strain 1
genoInP2 <- apply(calldata[,indPop2], 1, function(x){return(as.numeric(na.omit(unique(x))))})       # Genotypes in strain 2

# Differences, we need to take care of the hetrozygotes individuals:
# 2 2 2 2 2 2 2 versus 0 0 0 0 0 0 0 0 0
# 2 2 2 2 2 2 2 versus 0 0 0 0 0 1 0 0 0
# 2 2 1 2 1 1 2 versus 0 0 0 0 0 0 0 0 0
different <- rep(FALSE, nrow(calldata))
for(x in 1:nrow(calldata)){
  if(length(genoInP1[[x]])==1) different[x] = !(genoInP1[[x]] %in% genoInP2[[x]])
  if(length(genoInP2[[x]])==1) different[x] = !(genoInP2[[x]] %in% genoInP1[[x]])
  if(x %% 10000 == 0)cat("Done",x,"/",length(mightBeDifferent),"\n")
}

calldata <- calldata[which(different),]                                                         
write.table(calldata[,c(indPop1,indPop2)], file="Analysis/DifferentialGenotypes.txt", sep="\t")     # Different genotypes between the two RI strains

arrayAnnotation <- read.table("Annotation/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", header=TRUE, colClasses="character", na.string="---")

annotSubset <- which(arrayAnnotation[,"Probe.Set.ID"] %in% rownames(calldata))                      # Reduce the annotation file to contain only the SNPs
arrayAnnotation <- arrayAnnotation[annotSubset,]                                                    # different between the two RI strains

affyRelations <- c("splice-site", "intron", "nonsense", "missense", "synon", "frameshift", "aa-indel", "CDS", "UTR-3", "UTR-5", "exon", "upstream", "downstream")
# Helper function to count how many times a pattern is in a different string
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

write.table(relationCount, "Analysis/SNPeffectsInChickens.txt", sep="\t")
