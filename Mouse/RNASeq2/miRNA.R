#
# Correlation and differential correlation analysis for reciprocal crosses
#


setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross/Analysis")
alldata <- read.table(file=paste0("RPKM_norm_log_stats.txt"), sep="\t", header=TRUE)
rownames(alldata) <- alldata[,"ensembl_gene_id"]

allMIRs <- alldata[which(grepl("mir", alldata[,"gene_name"])),"gene_name"]


setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross/input")
miRNAs <- read.table("miRNA_Tph2.txt")

mirs <- gsub("R", "r", gsub("-5p","",gsub("-3p","",as.character(miRNAs[,1]))))  # Transform
ii <- which(as.character(alldata[,"gene_name"]) %in% mirs)  # Match

mirs <- gsub("mmu-mir-", "Mir", gsub("R", "r", gsub("-5p","",gsub("-3p","",as.character(miRNAs[,1])))))  # Transform
jj <- which(as.character(alldata[,"gene_name"]) %in% mirs)  # Match

write.table(alldata[c(ii,jj), ], "MIR_RPKM.txt", sep="\t")  # Collect