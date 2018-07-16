library("haplo.stats")

setwd("D:/Edrive/Goat/DNA/Sequencing/")
sampleIDs <- read.table("samples_description.txt", colClasses="character", sep="\t", header=FALSE,na.strings=c(""),row.names=1, check.names=FALSE)

setwd("D:/Edrive/Goat/DNA/Sequencing/SNPs")
snps <- read.table("SNPsInTranscriptsMAFGeno.txt", colClasses="character", sep="\t", header=TRUE,na.strings=c(""),check.names=FALSE)

snps <- t(snps)
snps <- snps[,which(snps["Synonymous",] == "Non")]

IDs <- unlist(lapply(strsplit(snps["ID",], ";"),"[",1))

chrs <- c("6", "8", "14")
names(chrs) <- c("ENA|CM004567|CM004567.1", "ENA|CM004569|CM004569.1", "ENA|CM004575|CM004575.1")

novelIDs <- paste0(chrs[snps["CHROM",which(is.na(IDs))]], "-", snps["POS",which(is.na(IDs))], snps["REF",which(is.na(IDs))], ">", snps["ALT",which(is.na(IDs))])

IDs[which(is.na(IDs))] <- novelIDs
colnames(snps) <- IDs

nonDuplicated = which(!duplicated(IDs))
snps <- snps[,nonDuplicated]

genes <- c("CSN1S1", "CSN2", "CSN1S2", "CSN3")

for(gene in genes){
  nSnps <- which(snps["inGene",] == gene)
  samples <- 6:38

  geneSNPs <- snps[samples, nSnps]
  geneSNPs[is.na(geneSNPs)] <- "XX"

  haploinput <- NULL
  for(x in 1:ncol(geneSNPs)){
    haploinput <- cbind(haploinput, matrix(unlist(strsplit(geneSNPs[,x], "")),length(samples),2,byrow=TRUE))
  }

  haploinput[haploinput == "X"] <- 0

  groups <- sampleIDs[rownames(snps[samples,]),1]
  groups[groups == "Wild"] <- "Ibex"

  haplogroups <- haplo.group(groups, haploinput, colnames(geneSNPs), haplo.em.control(max.iter = 100, min.posterior = 1e-06, verbose=TRUE))
  haplogroups = haplogroups[[1]]
  colnames(haplogroups) <- gsub("groups=", "", colnames(haplogroups))
  rownames(haplogroups) <- paste("Haplotype", LETTERS[1:nrow(haplogroups)])

  write.table(haplogroups, file=paste0("haplo_",gene,".txt"), sep="\t", na="")
}
