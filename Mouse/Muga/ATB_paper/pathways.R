setwd("D:/Edrive/Mouse/DNA/MegaMuga/Analysis")

pat_bfmi <- read.csv("PAT_BFMIgenesInfo.txt",sep="\t", header=TRUE, colClasses=c("character", rep("numeric",3), "character", "character"))
mat_bfmi <- read.csv("MAT_BFMIgenesInfo.txt",sep="\t", header=TRUE, colClasses=c("character", rep("numeric",3), "character", "character"))
pat_b6 <- read.csv("PAT_B6NgenesInfo.txt",sep="\t", header=TRUE, colClasses=c("character", rep("numeric",3), "character", "character"))
mat_b6 <- read.csv("MAT_B6NgenesInfo.txt",sep="\t", header=TRUE, colClasses=c("character", rep("numeric",3), "character", "character"))

pat_bfmi <- cbind(pat_bfmi, O = "PAT_BFMI")
mat_bfmi <- cbind(mat_bfmi, O = "MAT_BFMI")
pat_b6 <- cbind(pat_b6, O = "PAT_B6N")
mat_b6 <- cbind(mat_b6, O = "MAT_B6N")

allgenes <- rbind(pat_bfmi,mat_bfmi,pat_b6,mat_b6)
allgenes <- allgenes[-grep("predicted", allgenes[,"mgi_description"]),]
allgenes <- allgenes[-grep("RIKEN cDNA", allgenes[,"mgi_description"]),]
allgenes <- allgenes[-which(allgenes[,"mgi_description"] == ""),]

length(which(allgenes[, "O"] == "PAT_BFMI"))
length(which(allgenes[, "O"] == "MAT_BFMI"))
length(which(allgenes[, "O"] == "PAT_B6N"))
length(which(allgenes[, "O"] == "MAT_B6N"))

pat_bfmi <- allgenes[which(allgenes[, "O"] == "PAT_BFMI"),]
pat_bfmi <- pat_bfmi[with(pat_bfmi, order(chromosome_name, start_position)),]

mat_bfmi <- allgenes[which(allgenes[, "O"] == "MAT_BFMI"),]
mat_bfmi <- mat_bfmi[with(mat_bfmi, order(chromosome_name, start_position)),]

pat_b6 <- allgenes[which(allgenes[, "O"] == "PAT_B6N"),]
pat_b6 <- pat_b6[with(pat_b6, order(chromosome_name, start_position)),]

mat_b6 <- allgenes[which(allgenes[, "O"] == "MAT_B6N"),]
mat_b6 <- mat_b6[with(mat_b6, order(chromosome_name, start_position)),]

write.table(pat_bfmi, file="PAT_BFMI_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(mat_bfmi, file="MAT_BFMI_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(pat_b6, file="PAT_B6N_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(mat_b6, file="MAT_B6N_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)

cat("Unique genes in region: ", length(unique(allgenes[,1])), "\n")
cat("Overlap between PAT & MAT:, ", dim(allgenes)[1] - length(unique(allgenes[,1])), "\n")

snpdata <- read.table("D:/Edrive/Mouse/DNA/Sequencing/BFMI/20140515_VEP_BFMI860mm10.txt", sep = "\t", header=TRUE)

SNPsInTRD <- snpdata[which(as.character(snpdata[,"Gene_ID"]) %in% allgenes[, "ensembl_gene_id"] & grepl("CODING_REGION", snpdata[,"Region"]) & snpdata[,"FuncClass"] == "NON_SYNONYMOUS_CODING"), c(1:5,8,12,13,14,15,16,19)]

cat("Unique SNPs: ", length(unique(SNPsInTRD[,"SNP_ID"])), "\n")
cat("Unique Genes: ", length(unique(SNPsInTRD[,"Gene_ID"])), "\n")

snpgenes <- allgenes[which(allgenes[,1] %in% SNPsInTRD[,"Gene_ID"]),]
snpgenes <- snpgenes[with(snpgenes, order(chromosome_name, start_position)),]

pat_bfmi_snp <- snpgenes[which(snpgenes[, "O"] == "PAT_BFMI"),]
mat_bfmi_snp <- snpgenes[which(snpgenes[, "O"] == "MAT_BFMI"),]
pat_b6_snp <- snpgenes[which(snpgenes[, "O"] == "PAT_B6N"),]
mat_b6_snp <- snpgenes[which(snpgenes[, "O"] == "MAT_B6N"),]

write.table(pat_bfmi_snp, file="PAT_BFMI_NSyn_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(mat_bfmi_snp, file="MAT_BFMI_NSyn_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(pat_b6_snp, file="PAT_B6N_NSyn_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(mat_b6_snp, file="MAT_B6N_NSyn_Filtered.txt", sep = "\t", row.names=FALSE, quote=FALSE)



setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
results <- read.table("ChiSqScores_Incompatible.txt", sep="\t", check.names=FALSE)
rr <- apply(results,1,as.numeric)
rownames(rr) <- colnames(rr)
results <- rr

# Generate an overview plot
nTests <- (ncol(results) * ncol(results)) / 2
LODscores <- -log10(pchisq(results, 1, lower.tail=FALSE))
threshold <- -log10(0.05 / nTests)

significant <-  names(which(apply(LODscores,1,max,na.rm=TRUE) > threshold))

allRegions <- allRegions[significant, ]
regionnames <- rownames(allRegions)
results <- results[significant,significant]


which(results[1,] > threshold)


# Load all the known genes with MGI symbols
completegenome <- read.csv("D:/Edrive/Mouse/DNA/MegaMuga/Additional/MGI_Gene_Model_Coord.rpt", sep='\t', header=TRUE, check.names=FALSE)
completegenome <- completegenome[-grep("predicted", completegenome[,"3. marker symbol"]),]
completegenome <- completegenome[-grep("RIKEN cDNA", completegenome[,"3. marker symbol"]),]

# Permutation of the amount of non-synonymous SNPs
permSNPs <- NULL
permGenes <- NULL
for(x in 1:50000) {
  rselection <- as.character(completegenome[sample(nrow(completegenome), length(unique(allgenes[,1]))),2])
  shortlist <- snpdata[which(as.character(snpdata[,"Gene_Name"]) %in% rselection & grepl("CODING_REGION", snpdata[,"Region"]) & snpdata[,"FuncClass"] == "NON_SYNONYMOUS_CODING"), c(1:5,8,12,13,14,15,16,19)]
  permSNPs <- c(permSNPs, nrow(shortlist[-which(duplicated(shortlist[,"SNP_ID"])), ]))
  permGenes <- c(permGenes, nrow(shortlist[-which(duplicated(shortlist[,"Gene_ID"])), ]))
  cat("Done", x, "\n")
}
cat(permSNPs, sep = "\n", file = "permutations_SNPs.txt")
cat(permGenes, sep = "\n", file = "permutations_Genes.txt")

snpgenes[1:10,]

