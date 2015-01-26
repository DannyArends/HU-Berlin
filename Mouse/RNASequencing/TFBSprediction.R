# TFBSprediction.R - predict TFBS in upstream region of mouse
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
#
# Create a figure for the RNA sequencing data (pre-processed by MDC)

library(biomaRt)                                                                                                      # Biomart - Used to get upstream of genes
library(seqinr)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
RPKM <- read.table("Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")

RPKM <- RPKM[which(RPKM[,"Mean.B6NxBFMI860.12.L"] > 1 | RPKM[,"Mean.BFMI860.12xB6N.L"] > 1),]
RPKM <- RPKM[which(RPKM[,"tTest_F1"] < 0.05),]

mart <- useMart("ensembl", "mmusculus_gene_ensembl")
upstreamseqs <- getBM(attributes=c('gene_flank','start_position','end_position','chromosome_name','strand','ensembl_gene_id'), filters=c('ensembl_gene_id','upstream_flank'), values=list(ENSG = RPKM[,"ensembl_gene_id"], Upstream=2000), mart=mart, checkFilters=FALSE)

write.table(upstreamseqs, "Analysis/UpstreamSequences0.05.txt", sep = "\t")

cat("", file = "Analysis/UpstreamSequences0.05.fasta")
for(x in 1:nrow(upstreamseqs)){
  cat(">", upstreamseqs[x,"ensembl_gene_id"],"\n", sep = "", file = "Analysis/UpstreamSequences0.05.fasta", append=TRUE)
  cat(upstreamseqs[x,"gene_flank"],"\n", sep = "", file = "Analysis/UpstreamSequences0.05.fasta", append=TRUE)
}

### Run clustalW2 to align the sequences
### E:\Mouse\RNA\Sequencing\Reciprocal Cross B6 BFMI by MPI\Analysis\UpstreamSequences0.05.fasta

sequences <- read.alignment("Analysis/UpstreamSequences0.05.aln")


### Check for motives in the upstream of the aligned sequences
library(rtfbs)

ms <- read.ms("Analysis/UpstreamSequences0.05.fasta")
tmp <- groupByGC.ms(ms, 1)
for (i in 1:length(tmp)){
  mm <- build.mm(tmp[[i]], 40)
  sreal <- score.ms(tmp[[i]], pwm, mm)
  ss <- simulate.ms(mm, sum(lengths.ms(ms)))
  ssim <- score.ms(ss, pwm, mm)
  fdrMap <- calc.fdr(tmp[[i]], sreal, ss, ssim)
  makeFdrPlot(fdrMap)
  output.sites(sreal, fdrScoreMap=fdrMap, fdrThreshold = 0.05)
}