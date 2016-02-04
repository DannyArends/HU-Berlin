setwd("E:/Mouse/DNA/Sequencing/BFMI")

if(!file.exists("20140515_VEP_BFMI860mm10.txt")){
  library(openxlsx)
  snpdata <- NULL
  for(x in 1:22){
    snpdata <- rbind(snpdata, read.xlsx("20140515_VEP_BFMI860mm10.xlsx", sheet = x))
  }
  write.table(snpdata, "20140515_VEP_BFMI860mm10.txt", sep = "\t", row.names=FALSE)
}else{
  snpdata <- read.table("20140515_VEP_BFMI860mm10.txt", sep = "\t", header=TRUE)
}

setwd("E:/Mouse/DNA/MegaMuga/")
allgenes <- NULL                                                                                  # Get all genes in each regions
allgenes <- read.table("Analysis/Mat_B6NgenesInfo.txt", sep="\t", header=TRUE)
allgenes <- rbind(allgenes, read.table("Analysis/Mat_BFMIgenesInfo.txt", sep="\t", header=TRUE))
allgenes <- rbind(allgenes, read.table("Analysis/Pat_B6NgenesInfo.txt", sep="\t", header=TRUE))
allgenes <- rbind(allgenes, read.table("Analysis/Pat_BFMIgenesInfo.txt", sep="\t", header=TRUE))

regionsMat <- read.table("regions_matp0.01.txt", sep="\t", header=TRUE)
regionsPat <- read.table("regions_patp0.01.txt", sep="\t", header=TRUE)
allRegions <- rbind(cbind(regionsMat, origin = "M"), cbind(regionsPat, origin = "P"))
allRegions <- allRegions[with(allRegions, order(Chr, Start)), ]

regionnames <- paste0(allRegions[,"Chr"],":", 
          round(allRegions[,"Start"] / 1000000, 0),"-",
          round(allRegions[,"Stop"] / 1000000, 0), " ", 
          allRegions[,"origin"])

rownames(allRegions) <- regionnames

scores <- read.table("ChiSqScores_Incompatible.txt", sep="\t", check.names=FALSE)
scoresnum <- apply(scores,1,as.numeric)
colnames(scoresnum) <- colnames(scores)
rownames(scoresnum) <- rownames(scores)
LODscores <- -log10(pchisq(scoresnum, 1, lower.tail=FALSE))

# Which region the gene is in
whichRegion <- function(allgenes, ensgid, verbose = TRUE){
  ii <- allgenes[which(allgenes[,1] == ensgid),]
  cat(ensgid, ": Found", nrow(ii), "matching genes\n")
  gChr <- ii[1,"chromosome_name"];gLoc <- ii[1,"start_position"]  # Chr:Location
  cat(ensgid, ": ", gChr, "-", gLoc, "\n")
  rownames(allRegions[which(allRegions[,"Chr"] == gChr & as.numeric(allRegions[,"Stop"]) > gLoc ),])
}

genenames <- unique(as.character(allgenes[,"mgi_symbol"]))    # All genes in the regions

# How many genes in the TRD regions with non-synonimous SNPs
shortlist <- snpdata[which(as.character(snpdata[,"Gene_Name"]) %in% genenames & grepl("CODING_REGION", snpdata[,"Region"]) & snpdata[,"FuncClass"] == "NON_SYNONYMOUS_CODING"), c(1:5,8,12,13,14,15,16,19)]
write.table(shortlist[-which(duplicated(shortlist[,"SNP_ID"])), ], "AllRegionsNonSynonimousGenes.txt", sep = "\t", row.names=FALSE)
nrow(shortlist[-which(duplicated(shortlist[,"SNP_ID"])), ])

# 5% threshold for the genetic incompatibilities
threshold <- -log10(0.05 / ((ncol(LODscores) * ncol(LODscores)) / 2))

allgenes <- allgenes[which(allgenes[,1] %in% shortlist[,"Gene_ID"]),]
allgenes <- allgenes[-which(duplicated(allgenes[,5])),]
allgenes <- cbind(allgenes, region = NA)

# Find which region each gene is in
for(x in 1:nrow(allgenes)){
  allgenes[x, "region"] <- whichRegion(allgenes, allgenes[x,"ensembl_gene_id"])[1]
}

# List possible interactions between genes showing non-synonymous SNPs
results <- NULL
for(x in 1:nrow(allgenes)){
  gene1 <- as.character(allgenes[x, "mgi_symbol"])
  r1 <- allgenes[x,"region"]
  if(!is.na(r1)){
    ii <- which(LODscores[r1, ] > threshold)
    cat(gene1, ii,"\n")
    for(i in ii){
      r2 <- colnames(LODscores)[i]
      cat(r1, r2, LODscores[r1,r2], "\n")
      jj <- which(allgenes[,"region"] == r2)
      for(j in jj){
        gene2 <- as.character(allgenes[j,"mgi_symbol"])
        if(allgenes[x,"chromosome_name"] < allgenes[j,"chromosome_name"]){
          results <- rbind(results, cbind(allgenes[x,], allgenes[j,], LODscores[r1,r2]))
        }
        #cat(gene1, r1, gene2, r2, LODscores[r1,r2], "\n")
      }
    }
  }
}
results <- results[with(results, order(chromosome_name, start_position)), ]
colnames(results) <- c("ensembl_gene_id", "chr", "start", "end", "mgi_symbol", "mgi_description", "region", "ensembl_gene_id", "chr", "start", "end", "mgi_symbol", "mgi_description", "region", "LOD")

write.table(results,"InteractionsNonSynGenes.txt",sep="\t", row.names=FALSE)

# Load all the known genes with MGI symbols
completegenome <- read.table("Additional/MGI_Gene_Model_Coord.rpt", sep='\t')
# Permutation of the amount of non-synonymous SNPs
perms <- NULL
for(x in 1:10000) {
  rselection <- as.character(completegenome[sample(nrow(completegenome),1158),3])
  shortlist <- snpdata[which(as.character(snpdata[,"Gene_Name"]) %in% rselection & grepl("CODING_REGION", snpdata[,"Region"]) & snpdata[,"FuncClass"] == "NON_SYNONYMOUS_CODING"), c(1:5,8,12,13,14,15,16,19)]
  perms <- c(perms, nrow(shortlist[-which(duplicated(shortlist[,"SNP_ID"])), ]))
}

