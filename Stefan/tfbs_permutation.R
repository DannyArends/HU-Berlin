library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.ensGene)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

transcripts       <- transcriptsBy(TxDb.Mmusculus.UCSC.mm10.ensGene, by="gene")
transcriptsOnChr  <- lapply(transcripts, function(x){ return(as.character(x@seqnames@values))})
chrs              <- unique(unlist(transcriptsOnChr))[-grep("_", unique(unlist(transcriptsOnChr)))]
transcripts       <- transcripts[which(transcriptsOnChr %in% chrs)]
upstreams         <- getPromoterSeq(transcripts, Mmusculus, upstream=2000, downstream=1000)

genesDE    <- read.table("tfbs/papergenes.txt", colClasses="character", sep="\t", header=TRUE)

upgenes    <- as.character(unlist(genesDE[genesDE[,"Ratio_F1"] > 1,"ensembl_gene_id"])); upgenes    <- upgenes[-which(!upgenes %in% names(upstreams))]
downgenes  <- as.character(unlist(genesDE[genesDE[,"Ratio_F1"] < 1,"ensembl_gene_id"])); downgenes  <- downgenes[-which(!downgenes %in% names(upstreams))]

TFBSmatches <- function(WeightMatrix, selgenes, upstreams){
  motifTFBS <- round( 100 * WeightMatrix)

  nmatches <- NULL
  for(x in 1:length(selgenes)){
    if(length(matchPWM(motifTFBS, unlist(upstreams[selgenes[x]])[[1]], "90%")) > 0) nmatches <- c(nmatches, selgenes[x])
  }
  #cat("Matches =",length(nmatches),"\n") #,", p-value Up/Down = ", pUp, "/", pDown,"\n")
  return(nmatches)
}

TFBSperm <- function(ingenes, TFBSlist, nperm = 1000) {
  mREAL <- NULL
  for(TFBS in TFBSlist){
    mREAL <- c(mREAL, TFBSmatches(TFBS, ingenes, upstreams))
  }
  nmatches <- length(unique(mREAL))

  nulldistribution <- NULL
  for(p in 1:nperm){
    rangenes <- names(upstreams)[sample(length(upstreams),length(ingenes))]
    mPERM <- NULL
    for(TFBS in TFBSlist){ mPERM <- c(mPERM, TFBSmatches(TFBS, rangenes, upstreams)) }
    nulldistribution <- c(nulldistribution, length(unique(mPERM)))
    cat(paste0("Permutation ", p, ": ", length(unique(mPERM))),"\n")
  }
  pUp   <- (1 - (length(which(nulldistribution <= nmatches)) / length(nulldistribution)))
  pDown <- (1 - (length(which(nulldistribution >= nmatches)) / length(nulldistribution)))
  cat(nmatches," ", pUp," ", pDown,"\n")
  c(nmatches, nulldistribution)
}

### Glucocorticoid weight matrices

weights.GR <- rbind(
c(0.24,0.21,0.21,0.24,0.24,0.37,0.13,0.39,0.32,0.26,0.21,0,0,0.11,0.10,0.05,0.24,0.32,0.20),      #A
c(0.15,0.21,0.11,0.21,0.24,0.24,0.58,0.21,0.34,0.08,0.26,0,0,0.02,0.24,0.79,0.13,0.21,0.24),      #C
c(0.24,0.37,0.47,0.45,0.18,0.10,0.16,0.16,0.10,0.08,0.39,0,1,0,0.16,0.13,0,0.15,0.32),            #G
c(0.37,0.21,0.21,0.10,0.34,0.29,0.13,0.24,0.24,0.58,0.14,1,0,0.87,0.50,0.03,0.63,0.32,0.24))      #T
rownames(weights.GR) <- c("A","C","G","T")

resUpGC   <- TFBSperm(upgenes,   list(weights.GR)) # Up regulated genes  p < 0.05
resDownGC <- TFBSperm(downgenes, list(weights.GR)) # Down regulated genes  p < 0.05

### Estrogen weight matrices

weights.ER <- rbind(
c(0.3,0.2,0.6,0.3,0.15,0.15,0.15,0.55,0.2,0.15,0.1,0,0.05,1,0,0.05,0.1,0.35,0.4),         #A
c(0.2,0.45,0.2,0.05,0.1,0.15,0.5,0.1,0.45,0.3,0.15,0.05,0,0,0.9,0.9,0.35,0.05,0.25),      #C
c(0.2,0.15,0.15,0.55,0.55,0.2,0.15,0.2,0.15,0.15,0.45,0,0.95,0,0.1,0.05,0.1,0.4,0.15),    #G
c(0.3,0.2,0.05,0.1,0.2,0.5,0.2,0.15,0.2,0.4,0.3,0.95,0,0,0,0,0.45,0.2,0.2))               #T
rownames(weights.ER) <- c("A","C","G","T")

weights.ERE <- rbind(
c(0,0,0,0,1,0.25,0.25,0.25,0,0,1,0,0),      #A
c(0,0,0,1,0,0.25,0.25,0.25,0,0,0,1,1),      #C
c(1,1,0,0,0,0.25,0.25,0.25,0,1,0,0,0),      #G
c(0,0,1,0,0,0.25,0.25,0.25,1,0,0,0,0))      #T
rownames(weights.ERE) <- c("A","C","G","T")

weights.ESR1 <- rbind(
c(0.11,0.11,0.78,0.22,0,0,0,0.67,0.11,0.22,0.33,0.11,0.11,0.56,0,0,0.22,0.33),            #A
c(0.56,0.56,0.11,0,0,0,0.78,0,0.78,0.56,0.22,0.11,0,0.11,0.89,1,0.45,0.45),               #C
c(0.11,0.11,0.11,0.78,1,0,0.22,0.22,0.11,0.11,0.45,0.11,0.89,0.33,0,0,0,0.11),            #G
c(0.22,0.22,0,0,0,1,0,0.11,0,0.11,0,0.67,0,0,0.11,0,0.33,0.11))                           #T
rownames(weights.ESR1) <- c("A","C","G","T")

JASPER <- query(MotifDb, "Mmusculus")
Esrra  <- "Mmusculus-UniPROBE-Esrra.UP00079";      Esrra2 <- "Mmusculus-jolma2013-Esrra-2"
Esrrb  <- "Mmusculus-JASPAR_CORE-Esrrb-MA0141.1";  Esrrb2 <- "Mmusculus-JASPAR_2014-Esrrb-MA0141.2"

resUpEstro   <- TFBSperm(upgenes,   list(weights.ER, weights.ERE, weights.ESR1, JASPER[Esrra][[1]], JASPER[Esrra2][[1]], JASPER[Esrrb][[1]], JASPER[Esrrb2][[1]])) # Up regulated genes  p < 0.05
resDownEstro <- TFBSperm(downgenes, list(weights.ER, weights.ERE, weights.ESR1, JASPER[Esrra][[1]], JASPER[Esrra2][[1]], JASPER[Esrrb][[1]], JASPER[Esrrb2][[1]])) # Down regulated genes  p < 0.05

### Androgen weight matrices
weights.ARE <- rbind(
c(1,0,1,1,0,1,0.25,0.25,0.25,0,0,0,1,0,0),  #A
c(0,0,0,0,1,0,0.25,0.25,0.25,1,0,0,0,1,0),  #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0),  #G
c(0,0,0,0,0,0,0.25,0.25,0.25,0,0,1,0,0,1))  #T
rownames(weights.ARE) <- c("A","C","G","T")

weights.ARE2 <- rbind(
c(1,0,0,1,0,0,0.25,0.25,0.25,1,0,1,1,0,0),  #A
c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,1),  #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0),  #G
c(0,0,1,0,0,1,0.25,0.25,0.25,0,0,0,0,0,0))  #T
rownames(weights.ARE2) <- c("A","C","G","T")

weights.ADR3 <- rbind(
c(1,0,1,1,0,1,0.25,0.25,0.25,1,0,1,1,0,1),  #A
c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,0),  #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,0),  #G
c(0,0,0,0,0,0,0.25,0.25,0.25,0,0,0,0,0,0))  #T
rownames(weights.ADR3) <- c("A","C","G","T")

weights.IR3 <- rbind(
c(1,0,1,1,0,1,0.25,0.25,0.25,0,0,0,0,0,0), #A
c(0,0,0,0,1,0,0.25,0.25,0.25,0,0,0,0,1,0), #C
c(0,1,0,0,0,0,0.25,0.25,0.25,0,1,0,0,0,1), #G
c(0,0,0,0,0,0,0.25,0.25,0.25,1,0,1,1,0,0)) #T
rownames(weights.IR3) <- c("A","C","G","T")

resUpAG   <- TFBSperm(upgenes,   list(weights.ARE, weights.ARE2, weights.ADR3, weights.IR3, query(MotifDb, "Mmusculus")["Mmusculus-jolma2013-Ar"][[1]])) # Up regulated genes  p < 0.05
resDownAG <- TFBSperm(downgenes, list(weights.ARE, weights.ARE2, weights.ADR3, weights.IR3, query(MotifDb, "Mmusculus")["Mmusculus-jolma2013-Ar"][[1]])) # Down regulated genes  p < 0.05
