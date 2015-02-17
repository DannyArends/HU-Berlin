
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
GTF <- read.table("GTF/Mus_musculus.GRCm38.76.gtf", sep="\t")                                                       # Gene models

geneids <- gsub("gene_id ", "", unlist(lapply(strsplit(as.character(GTF[,9]),";"),"[",1)))
bbs7 <- GTF[which(geneids == "ENSMUSG00000037325"),]    # Only BBS7

cat(readLines("GTF/Mus_musculus.GRCm38.76.gtf",n=5),file="GTF/Bbs7.gtf",sep="\n")
write.table(bbs7,file="GTF/Bbs7.gtf", append=TRUE, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

## Create exon level expression dataset from the RNA-Seq data

files <- c("4868_GCCAAT_L001_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5067_ATCACG_L004_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5068_GAGTGG_L004_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5069_AGTCAA_L004_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5070_CGATGT_L005_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5071_CCGTCC_L005_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5072_TAGCTT_L005_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5073_TTAGGC_L006_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5074_GATCAG_L006_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5075_ATGTCA_L006_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5076_TGACCA_L007_P_trimmed.aligned.sorted.realigned.recalibrated.bam",
"5077_ACTTGA_L007_P_trimmed.aligned.sorted.realigned.recalibrated.bam")


for(file in files){
  cat("nohup htseq-count -r name -f bam -s no -a 0 -t exon -i exon_id -q Analysis/", file," GTF/Mus_musculus.GRCm38.76.gtf > Analysis/",substr(file,0,4),".exons &\n",sep="")
}


library("preprocessCore")

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
alldata <- 0
for(mfile in dir("ExonLevel")){
  fdata <- read.table(paste0("ExonLevel/",mfile),sep="\t",skip=2)
  colnames(fdata) <- paste0(c("name","reads"),"_", substr(mfile,1,4))
  alldata <- cbind(alldata, fdata)
}

rawreads <- alldata[, grep("reads", colnames(alldata))]
rownames(rawreads) <- alldata[,2]

rawreads[rawreads == 0] <- NA                                                              # Change 0 reads to NA, so we can do quantile normalisation
rawreadsLog2 <- log2(rawreads)
rawreadsQnorm <- normalize.quantiles(as.matrix(rawreadsLog2))
rawreadsQnorm[is.na(rawreadsQnorm)] <- 0                                                    # Change back the NAs to 0 reads

boxplot(rawreadsQnorm)

