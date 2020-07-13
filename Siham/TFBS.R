#
# TFBS analysis Siham
#
library(MotifDb)
library(seqLogo)

setwd("D:/Edrive/Goat/DNA/Sequencing")

chromosomes <- c(6,8,14)
names(chromosomes) <- c("ENA|CM004567|CM004567.1", "ENA|CM004569|CM004569.1", "ENA|CM004575|CM004575.1")

genes <- read.table("genes_with_sequence.txt", sep = "\t")
genes <- genes[1:4,]
snps <- read.table("SNPs/filteredSNPs.txt", sep = "\t", header=TRUE)
snps <- cbind(Chr = chromosomes[snps[, "CHROM"]], snps)

#genome <- readLines("Genome/CM004567.fasta.txt")[-1]
#wholegenome <- paste0(genome,collapse="")
#wholegenomeL <- strsplit(wholegenome, "")[[1]]
#genes <- cbind(genes, gSeq = NA)
#genes[1, "gSeq"] = paste0(wholegenomeL[genes[1,4]:genes[1,5]], collapse="")
#genes[2, "gSeq"] = paste0(wholegenomeL[genes[2,4]:genes[2,5]], collapse="")
#genes[3, "gSeq"] = paste0(wholegenomeL[genes[3,4]:genes[3,5]], collapse="")
#genes[4, "gSeq"] = paste0(wholegenomeL[genes[4,4]:genes[4,5]], collapse="")

mouseTFBS <- query(MotifDb, "Mmusculus")

res <- vector("list", length(mouseTFBS))
names(res) <- names(mouseTFBS)
for(x in 1:nrow(genes)){
  tfbslist <- c()
  for(tfbs in 1:length(mouseTFBS)){
    #seqLogo(mouseTFBS[[tfbs]])
    TFBSjaspar <- round(100 * mouseTFBS[[tfbs]])
    hitsref <- matchPWM(TFBSjaspar, as.character(genes[x, 6]), "99%")
    if (length(hitsref) > 0) {
      for (y in 1:length(hitsref)) {
        pS <- start(hitsref[y]) + genes[x,4]
        pE <- start(hitsref[y]) + width(hitsref[y]) + genes[x,4]
        ii <- which(snps[, "Chr"] == genes[x,"chr"] & snps[, "POS"] >= pS & snps[, "POS"] <= pE)
        if(length(ii) > 0){
          cat("SNP", ii, "in TFBS", tfbs, "\n")
          for(i in ii){
            snpPos <- (snps[i, "POS"] - genes[x,4]) + 1
            bp <- strsplit(as.character(genes[x, 6]), "")[[1]][snpPos]
            if(snps[i, "REF"] != bp) stop(paste0(snps[i, "REF"], "!=", bp, "\n"))
            
          }
        }
      }
    }
  }
}

