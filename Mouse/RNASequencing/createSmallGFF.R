
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")
GTF <- read.table("GTF/Mus_musculus.GRCm38.76.gtf", sep="\t")                                                       # Gene models
GTF <- GTF[GTF[,3] == "exon",]

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
colnames(rawreads) <- gsub("reads_","",colnames(rawreads))

rawreads <- rawreads[apply(rawreads, 1, sum) > 100, ]

rawreads[rawreads == 0] <- NA                                                              # Change 0 reads to NA, so we can do quantile normalisation
rawreadsQnorm <- normalize.quantiles(as.matrix(rawreads))
rawreadsQnorm[is.na(rawreadsQnorm)] <- 0                                                    # Change back the NAs to 0 reads
rownames(rawreadsQnorm) <- rownames(rawreads)
colnames(rawreadsQnorm) <- gsub("reads_","",colnames(rawreads))

boxplot(rawreadsQnorm)

BFMI  <- c("4868", "5067")
B6N   <- c("5068","5069")
mB6N  <- c("5070","5071","5072")
mBFMI <- c("5073","5074","5075")

pvalues_P <- apply(rawreadsQnorm[,c(BFMI,B6N)], 1, function(x){ res <- NA; tryCatch(res <- t.test(x[BFMI], x[B6N])$p.value, error = function(e) e); return(res) })
pvalues_M <- apply(rawreadsQnorm[,c(mBFMI,mB6N)], 1, function(x){ res <- NA; tryCatch(res <- t.test(x[mBFMI], x[mB6N])$p.value, error = function(e) e); return(res) })

mean_mBFMI <- apply(rawreadsQnorm[,c(mBFMI,mB6N)], 1, function(x){ mean(x[mBFMI]) })
mean_mB6N  <- apply(rawreadsQnorm[,c(mBFMI,mB6N)], 1, function(x){ mean(x[mB6N]) })
ratios_M <- log2(apply(rawreadsQnorm[,c(mBFMI,mB6N)], 1, function(x){ mean(x[mBFMI]) / mean(x[mB6N]) }))
FoldChange <- apply(rawreadsQnorm[,c(mBFMI,mB6N)], 1, function(x){ mean(x[mBFMI]) / mean(x[mB6N]) })

expressed <- mean_mBFMI > 50 | mean_mB6N > 50
rawreadsdiff <- cbind(rawreadsQnorm, FoldChange, Log2FC = ratios_M, Pvalue = pvalues_M)[which(expressed & abs(ratios_M) > 1.5),]
rawreadsdiff

#rawreadsdiff <- cbind(rawreadsQnorm,pvalues_M)[which(expressed & abs(pvalues_M) < 0.0001),]
#rawreadsdiff

annotation <- NULL
for(x in rownames(rawreadsdiff)){
  mindex <- grep(x,as.character(GTF[,9]))[1]
  msplit <- unlist(strsplit(as.character(GTF[mindex, 9]),";"))
  geneID <- gsub("gene_id ","", msplit[1])
  exonNUM <- gsub(" exon_number ","", msplit[3])
  geneName <- gsub(" gene_name ","", msplit[4])
  annotation <- rbind(annotation, c(x, geneID, geneName, exonNUM, as.character(GTF[mindex, 1]),GTF[mindex, 4],GTF[mindex, 5]))
}
colnames(annotation) <- c("exon_id","gene_id", "mgi_symbol","exon_num","chromosome", "exon_start", "exon_stop")
cat(unique(annotation[,"gene_id"]),sep="\n") # For easy GO / Over representation analysis

write.table(cbind(annotation, rawreadsdiff), file="Exon_differences_ratio_above_1.5.txt", sep = "\t", quote = FALSE)


