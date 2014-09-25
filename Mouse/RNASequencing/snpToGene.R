# snpToGene.R - Analyze the SNPs and indels called by the GenomeAnalysisToolKit, and summarize them into genes
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Sep, 2014

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
matB6Nsnps  <- read.csv("maternalB6snps_5reads.txt", sep="\t", header=TRUE)                                   # SNPs detected in the maternal B6N F1 cross
matBFMIsnps <- read.csv("maternalBFMIsnps_5reads.txt", sep="\t", header=TRUE)                                 # SNPs detected in the maternal BFMI F1 cross
RPKM        <- read.csv("Analysis/BFMI_RPKM_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")   # RPKM values from RNA-Seq

GTF <- read.table("GTF/Mus_musculus.GRCm38.76.gtf", sep="\t")                                                   # Gene models
EXONS <- GTF[which(GTF[,3]=="exon"),]

datastr <- strsplit(as.character(EXONS[,9]), "; ")
lengths <- unlist(lapply(datastr, length))

addInfo <- matrix(NA, length(datastr), 2)
for(x in 1:length(datastr)){ addInfo[x, ] <- c(strsplit(datastr[[x]][1]," ")[[1]][2], strsplit(datastr[[x]][lengths[x]]," ")[[1]][2]); }

EXONS <- cbind(addInfo, EXONS)
uniqueGenes <- unique(EXONS[,1])

geneExonsCoupling <- vector("list", length(uniqueGenes))                                                        # Create a list which enumerate the exons per gene
x <- 1
for(gene in uniqueGenes){
  geneExons <- which(EXONS[,1] == gene)
  possibleExons <- EXONS[geneExons[!duplicated(EXONS[geneExons,2])], -11]
  geneExonsCoupling[[x]] <- possibleExons
  cat(x, "Gene", gene, "has", nrow(geneExonsCoupling[[x]]), "/", length(geneExons),"Exons\n")
  x <- x + 1
}

summarizeImprintedSNPsInGene <- function(uniqueGenes, geneExonsCoupling, direction){
  geneSNPCoupling <- vector("list", length(uniqueGenes))
  x <- 1
  for(gene in geneExonsCoupling){                                                                               # Summarize all SNPs from a direction per gene
    onChr <- as.character(direction[,"Chr"]) == as.character(unique(gene[,3]))
    for(exon in 1:nrow(gene)){
      exonStart <- as.numeric(gene[exon,6])
      exonEnd <- as.numeric(gene[exon,7])
      inEXON <- which(direction[onChr,"Loc"] >= exonStart & direction[onChr,"Loc"] <= exonEnd)
      if(length(inEXON) > 0){
        geneSNPCoupling[[x]] <- rbind(geneSNPCoupling[[x]], cbind(direction[onChr,][inEXON,], Exon=exon))
        geneSNPCoupling[[x]] <- geneSNPCoupling[[x]][!duplicated(geneSNPCoupling[[x]][,"ID"]),]
      }
    }
    cat(x, "found", nrow(geneSNPCoupling[[x]]), "SNPs in gene\n")
    x <- x + 1
  }
  return(geneSNPCoupling)
  #mmatrix <- matrix(NA, length(geneSNPCoupling), 3)                                                             # Create an output matrix with 3 columns
  #for(x in 1:length(geneSNPCoupling)){ 
  #  mmatrix[x,1] <- as.character(uniqueGenes[x])                                                                # Ensembl geneID
  #  mmatrix[x,2] <- mean(geneSNPCoupling[[x]][,"ImprintingScore"]);                                             # Mean imprinting score across the gene
  #  origin <- names(which.max(table(geneSNPCoupling[[x]][,"Origin"])))
  #  if(!is.null(origin)){ mmatrix[x,3] <- origin; }                                                             # Take the one which occurs most as the origin
  #}
  #colnames(mmatrix) <- c("geneID", "ImprintingScore", "Origin")
  #return(mmatrix)
}

BFMIsummary <- summarizeImprintedSNPsInGene(uniqueGenes, geneExonsCoupling, matBFMIsnps)
names(BFMIsummary) <- uniqueGenes
B6Nsummary  <- summarizeImprintedSNPsInGene(uniqueGenes, geneExonsCoupling, matB6Nsnps)
names(B6Nsummary) <- uniqueGenes

imputeReferenceASE <- function(BFMIsummary, B6Nsummary, RPKM, RPKMcutoff = 3, ASEcutoff = 0.35){
  for(x in 1:length(BFMIsummary)){
    if(!is.null(BFMIsummary[[x]]) && is.null(B6Nsummary[[x]])){
      if(mean(BFMIsummary[[x]][,"ImprintingScore"]) >= ASEcutoff && RPKM[which(RPKM[,"ensembl_gene_id"] == names(BFMIsummary[x])),"Mean.B6NxBFMI860.12.L"] >= RPKMcutoff){
        B6Nsummary[[x]] <- BFMIsummary[[x]]
        for(snp in 1:nrow(B6Nsummary[[x]])){
          for(origin in c("Origin1", "Origin2", "Origin3")){
            if(B6Nsummary[[x]][snp,origin]=="BFMI"){ B6Nsummary[[x]][snp,origin] <- "B6N"; }else{ B6Nsummary[[x]][snp,origin] <- "BFMI"; }
          }
          B6Nsummary[[x]][snp,"ImprintingScore"] <- 1
          for(column in c("R1", "A1", "R2", "A2", "R3", "A3")){ B6Nsummary[[x]][snp, column] <- "?"; }
        }
        cat(x,"imputed\n");
      }
    }
    if(is.null(BFMIsummary[[x]]) && !is.null(B6Nsummary[[x]])){
      if(mean(B6Nsummary[[x]][,"ImprintingScore"]) >= ASEcutoff && RPKM[which(RPKM[,"ensembl_gene_id"] == names(BFMIsummary[x])),"Mean.BFMI860.12xB6N.L"] >= RPKMcutoff){
        BFMIsummary[[x]] <- B6Nsummary[[x]]
        for(snp in 1:nrow(BFMIsummary[[x]])){
          for(origin in c("Origin1", "Origin2", "Origin3")){
            if(BFMIsummary[[x]][snp,origin]=="BFMI"){ BFMIsummary[[x]][snp,origin] <- "B6N"; }else{ BFMIsummary[[x]][snp,origin] <- "BFMI"; }
          }
          BFMIsummary[[x]][snp,"ImprintingScore"] <- 1
          for(column in c("R1", "A1", "R2", "A2", "R3", "A3")){ BFMIsummary[[x]][snp, column] <- "?"; }
        }
        cat(x,"imputed\n");
      }
    }
  }
  return(list(BFMIsummary, B6Nsummary))
}


getShortList <- function(CROSSsummary, cutoff = 0.35){
  v <- NULL
  for(x in 1:length(CROSSsummary)){
    v <- c(v, mean(CROSSsummary[[x]][,"ImprintingScore"]))
  }
  hist(v)
  return(CROSSsummary[which(v > cutoff)])
}

BFMIase <- getShortList(BFMIsummary)
B6Nase <- getShortList(B6Nsummary)

BFMIsummary <- BFMIsummary[-which(is.na(BFMIsummary[,2])),]
B6Nsummary  <- B6Nsummary[-which(is.na(B6Nsummary[,2])),]

ordering <- match(BFMIsummary[,"geneID"], RPKM[,"ensembl_gene_id"])
write.table(cbind(RPKM[ordering,], BFMIsummary), "Imprinted_matBFMIsnps_5reads.txt", sep="\t", quote=FALSE, row.names=FALSE)

ordering <- match(B6Nsummary[,"geneID"], RPKM[,"ensembl_gene_id"])
write.table(cbind(RPKM[ordering,], B6Nsummary), "Imprinted_matB6Nsnps_5reads.txt", sep="\t", quote=FALSE, row.names=FALSE)
