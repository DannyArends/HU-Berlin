# investigate missing homozygous loci
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

library(biomaRt)

setwd("D:/Edrive/Mouse/DNA/MegaMuga/inputF2")
genotypes <- read.table(file="cleaned_genotypes_F2_All.txt", sep='\t')
founder.geno <- read.table(file="cleaned_genotypes_founders.txt", sep='\t')
map <- read.table(file="cleaned_map_25MbGap.txt", sep='\t')
phenotypes <- read.table(file="cleaned_phenotypes_F2.txt", sep='\t')

# Chromosome variables
chrs <- 1:20
names(chrs) <- c(1:19, "X")

mmu04910 <- read.table("mmu04910_Insulin signaling pathway.txt",sep='\t', row.names=NULL, header = FALSE, fill = TRUE)
mmu04068 <- read.table("mmu04068_FoxO signaling pathway.txt",sep='\t', row.names=NULL, header = FALSE, fill = TRUE)
mmu04922 <- read.table("mmu04922_Glucagon signaling pathway.txt",sep='\t', row.names=NULL, header = FALSE, fill = TRUE)

allgenes <- readLines("Mus_musculus.GRCm38.91.gff3")
allgenes <- allgenes[which(grepl("\tgene\t", allgenes))]
agm <- matrix(unlist(lapply(allgenes, strsplit, "\t")), length(allgenes), 9, byrow=TRUE)
geneids = unlist(lapply(agm[,9], function(x){ gsub("ID=gene:(.*?);(.+)", "\\1", x)}))
genesymbols = unlist(lapply(agm[,9], function(x){ gsub("(.+)Name=(.*?);(.+)", "\\2", x)}))
genedescriptions = unlist(lapply(agm[,9], function(x){ gsub("(.+)description=(.*?);(.+)", "\\2", x)}))
genedescriptions = unlist(lapply(strsplit(genedescriptions, " [Source", fixed=TRUE),"[",1))
genedescriptions <- gsub("%2C", ",", genedescriptions)
genedescriptions <- gsub("%3B", ";", genedescriptions)
genedescriptions <- gsub("%26", "&", genedescriptions)
agm <- cbind(ensembl_gene_id = geneids, gene_symbol = genesymbols, gene_description = genedescriptions, chromosome_name = agm[,1], start_position = agm[,4], end_position = agm[,5], strand = agm[,7])
agm <- agm[-which(grepl("predicted gene", agm[, "gene_description"])),]
agm <- agm[-which(grepl("RIKEN cDNA", agm[, "gene_description"])),]
agm <- agm[which(agm[,"chromosome_name"] %in% names(chrs)),]
agm <- agm[-which(duplicated(agm[,"gene_symbol"])),]
rownames(agm) <- agm[,"gene_symbol"]

kegg <- rbind(mmu04910, mmu04068, mmu04922)
kegg <- kegg[!duplicated(kegg[,2]),]
rownames(kegg) <- kegg[,2]
kegg <- kegg[,-c(2, 5)]

bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("external_gene_name"), values = rownames(kegg), mart = bio.mart)

HWE <- function(obs, expected = NULL){
  obs[is.na(obs)] <- 0
  if(is.null(expected)){
    P <- (2 * obs[1] + obs[2]) / (2 * ( obs[1] + obs[2] + obs[3]))
    Q <- (1 - P)
    n <- sum(obs)
    expAA <- (P^2) * n; expH <- 2 * P * Q * n; expBB <- (Q^2) * n
    expected <- c(expAA, expH, expBB)
  }
  chisq <- ((obs[1] - expected[1])^2 / expected[1]) + ((obs[2] - expected[2])^2 / expected[2]) + ((obs[3] - expected[3])^2 / expected[3])
  pval <- pchisq(chisq, 1, lower.tail = FALSE)
  return(pval)
}

res.biomart <- res.biomart[which(res.biomart[,"chromosome_name"] %in% names(chrs)),]
rownames(res.biomart) <- res.biomart[,"external_gene_name"]

kegg <- kegg[which(rownames(kegg) %in% rownames(res.biomart)),]
kegg <- cbind(kegg, res.biomart[rownames(kegg),])
colnames(kegg)[1:3] <- c("kegg_id", "kegg_description", "ko_id")

# Group sizes
gttables <- apply(genotypes, 1, table)

nA <- unlist(lapply(gttables,"[", "A"))
nH <- unlist(lapply(gttables,"[", "H"))
nB <- unlist(lapply(gttables,"[", "B"))
names(nA) <- names(nH) <- names(nB) <- rownames(genotypes)

#kegg <- agm

resultsFULL <- vector("list", nrow(kegg))
names(resultsFULL) <- rownames(kegg)
results <- NULL
for(x in 1:nrow(kegg)){
  moffset <- 0
  region <- c()
  while(length(region) < 4){
    region <- rownames(map[map[,"Chr"] == kegg[x, "chromosome_name"] & 
                      as.numeric(map[,"Mb_NCBI38"]) > (as.numeric(kegg[x, "start_position"]) - moffset) & 
                      as.numeric(map[,"Mb_NCBI38"]) < (as.numeric(kegg[x, "end_position"]) + moffset),])
    moffset <- moffset + 50000
  }
  regiongts <- cbind(map[region,1:2], nA = nA[region], nH = nH[region], nB = nB[region])
  foundergts <- founder.geno[region,]
  regiongts.founders <- regiongts
  toSwap <- which(foundergts[,1] == "B")
  if(length(toSwap) > 0){
    temp <- regiongts.founders[toSwap, "nB"]
    regiongts.founders[toSwap, "nB"] <- regiongts.founders[toSwap, "nA"]
    regiongts.founders[toSwap, "nA"] <- temp
  }
  colnames(regiongts.founders) <- c(colnames(map)[1:2], "BFMI", "H", "B6N")
  resultsFULL[[x]] <- list(moffset = as.character(moffset), size = nrow(regiongts), regiongts = regiongts, foundergts = foundergts, regiongts.founders = regiongts.founders)
  rGT <- round(apply(regiongts.founders[,c("BFMI","H", "B6N")], 2, function(x){ sum(as.numeric(x), na.rm=TRUE)})/nrow(regiongts.founders),0)
  status = "S"
  if(rGT[1] < 25) status = "B6N"
  if(rGT[3] < 25) status = "BFMI"
  results <- rbind(results, c(snpOffset = format(moffset, scientific=F), nSNP = nrow(regiongts), rGT, Status = status))
  if(x %% 100 == 0) cat("Done ",x,"/",nrow(kegg),"\n")
}
rownames(kegg) <- kegg[, "gene_symbol"]

#rownames(results) <- rownames(kegg)
results <- data.frame(results)
results <- cbind(kegg[rownames(results),c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")], results)
results$chromosome_name <- factor(results$chromosome_name, levels=names(chrs))

affected <- results[which(results[,"Status"] != "S"),]
affectedAll <- affected[with(affected, order(Status, chromosome_name, start_position)), ]


plot(c(min(map[onChr3, "Mb_NCBI38"]), max(map[onChr3, "Mb_NCBI38"])), c(0, 380), t = 'n')
points(map[onChr3, "Mb_NCBI38"], nA[onChr3], t = 'l', col=2)
points(map[onChr3, "Mb_NCBI38"], nH[onChr3], t = 'l', col=1, lwd=2)
points(map[onChr3, "Mb_NCBI38"], nB[onChr3], t = 'l', col=3)

ps <- apply(genotypes, 1, function(x){
  obsA <- length(which(as.character(x) == "A"))
  obsH <- length(which(as.character(x) == "H"))
  obsB <- length(which(as.character(x) == "B"))
  pval <- HWE(c(obsA, obsH, obsB), c(87,158,98))
})


plot(c(0, length(onChr3)), c(0, 380), t = 'n')
points(nA[onChr3], t = 'l', col=2)
points(nH[onChr3], t = 'l', col=1, lwd=2)
points(nB[onChr3], t = 'l', col=3)