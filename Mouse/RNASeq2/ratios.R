#
# Get the ratio's and test for difference
#

# After SNP calling

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

setwd("/home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/Analysis")
vcfdata         <- read.table("population.vcf", header = TRUE, colClasses="character")
write.table(cbind(vcfdata[,"CHROM"],vcfdata[,"POS"]),"SNPlocations.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

bamfiles        <- list.files(pattern=".aligned.sorted.realigned.bam$")
reference       <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.dna"
reference.fa    <- paste0(reference, ".fa")        #!

for(x in bamfiles){
  cat(x,"\n")
  execute(paste0("/home/neubert/Keane/samtools-1.2/samtools mpileup -uv -t DV -t DP -l SNPlocations.txt -f ", reference.fa, " ", x, " | bcftools call -c - > ",paste0(x,".vcf")),paste0(x,".vcf"))
}

dp4datafull <- NA
for(x in bamfiles){
  cat(x," - ")
  vcfdata <- read.table(paste0(x, ".vcf"),sep="\t",colClasses="character")
  vcfdata <- vcfdata[-which(grepl("INDEL", vcfdata[,8])),]
  write.table(cbind(vcfdata[, 1], vcfdata[, 2], vcfdata[, 5],"\t", gsub(",","\t", unlist(lapply(vcfdata[, 8], function(x){ sub(".*?DP4=(.*?);.*", "\\1", x)} )))), paste0(x,".dp4"), row.names=FALSE, col.names=FALSE, quote=FALSE)

  dp4data <- read.table(paste0(x,".dp4"),sep="\t")
  dp4data <- cbind(dp4data, dp4data[,2] + dp4data[,3])
  colnames(dp4data)[6] <- paste0(substr(x,1,4),"_Ref")
  dp4data <- cbind(dp4data, dp4data[,4] + dp4data[,5])
  colnames(dp4data)[7] <- paste0(substr(x,1,4),"_Alt")
  dp4data <- cbind(dp4data, Total = dp4data[,paste0(substr(x,1,4),"_Ref")] + dp4data[, paste0(substr(x,1,4),"_Alt")])
  dp4data <- dp4data[dp4data[,"Total"] > 5,]
  if(is.na(dp4datafull)){
    dp4datafull <- dp4data[,c(1,6,7)]
  }else{
    dp4datafull <- dp4datafull[which(dp4datafull[,1] %in% dp4data[,1]),]
    dp4data <- dp4data[which(dp4data[,1] %in% dp4datafull[,1]),]
    dp4datafull <- cbind(dp4datafull, dp4data[,c(6,7)])
  }
  cat(nrow(dp4datafull),"\n")
}
write.table(dp4datafull, "allsamples.dp4",sep="\t",quote=FALSE, row.names=FALSE)

#setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
reads <- read.table("allsamples.dp4",sep="\t",header=TRUE)

c("5070","5071","5072")
c("5073","5074","5075")

