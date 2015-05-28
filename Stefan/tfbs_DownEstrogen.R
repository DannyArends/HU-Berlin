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

downgenes  <- as.character(unlist(genesDE[genesDE[,"Ratio_F1"] < 1,"ensembl_gene_id"])); downgenes  <- downgenes[-which(!downgenes %in% names(upstreams))]

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
WeightMatrices <- list(weights.ER, weights.ERE, weights.ESR1, JASPER[Esrra][[1]], JASPER[Esrra2][[1]], JASPER[Esrrb][[1]], JASPER[Esrrb2][[1]])

plot(c(1,3000),c(1,length(downgenes)), t = 'n')
for(wm in 1:length(WeightMatrices)){
  motifTFBS <- round( 100 * WeightMatrices[[wm]])
  for(x in 1:length(downgenes)){
    matches <- matchPWM(motifTFBS, unlist(upstreams[downgenes[x]])[[1]], "90%")
    if(length(matches) > 0){
      for(y in 1:length(matches)){
        points(c(start(matches[y]),end(matches[y])), c(x,x), t = 'o', lwd = 2, col=wm+3)
      }
    }
  }
}

