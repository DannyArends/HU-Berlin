#
# Mpileup of SNP data obtained after using BCFtools to call SNPs in the whole population
#
#

# After SNP calling
setwd("/home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/Analysis")
matB6N    <- c("5070","5071","5072") ; matBFMI   <- c("5073","5074","5075") ; BFMI      <- c("4868", "5067") ; B6N       <- c("5068", "5069")

samtools.exec       <- "samtools"
reference       <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.dna"
reference.fa    <- paste0(reference, ".fa")        #!

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

# Load in the VCF data created by BCFtools
vcfdata         <- read.table("population.vcf", header = TRUE, colClasses="character")
write.table(cbind(vcfdata[,"CHROM"],vcfdata[,"POS"]),"SNPlocations.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

# Samtools mpileup of each bamfile at the locations of the SNPs detected, to obtain DP4 measurements
bamfiles        <- list.files(pattern = ".aligned.sorted.realigned.recalibrated.bam$")
for(x in bamfiles){
  indexed.bamfile <-  gsub(".bam", ".bai", x)                           # index the bamfile
  execute(paste0(samtools.exec, " index ", x), indexed.bamfile)         # index the bamfile
  cat("Mpileup of file:", x, "\n")
  execute(paste0("/home/neubert/Keane/samtools-1.2/samtools mpileup -uv -t DV -t DP -l SNPlocations.txt -f ", reference.fa, " ", x, " | bcftools call -c - > ",paste0(x,".vcf")),paste0(x,".vcf"))
}

getDP4 <- function(x){ 
  sub(".*?DP4=(.*?);.*", "\\1", x)
}

# Load in the created VCF files, and extract the DP4 for all locations which exist across all files
vcffiles        <- list.files(pattern = ".aligned.sorted.realigned.bam.vcf$")
dp4datafull <- NA
for(x in vcffiles){
  cat(x," - ")
  vcfdata <- read.table(x, sep="\t",colClasses="character")       # Read the VCF data
  vcfdata <- vcfdata[-which(grepl("INDEL", vcfdata[,8])),]        # Remove the indels
  write.table(cbind(vcfdata[, 1], vcfdata[, 2], vcfdata[, 4], "\t", gsub(",", "\t", unlist(lapply(vcfdata[, 8], getDP4)))), paste0(x,".dp4"), row.names=FALSE, col.names=FALSE, quote=FALSE)

  dp4data <- read.table(paste0(x,".dp4"),sep="\t")
  dp4data <- cbind(dp4data, dp4data[,2] + dp4data[,3])
  colnames(dp4data)[6] <- paste0(substr(x,1,4),"_Ref")
  dp4data <- cbind(dp4data, dp4data[,4] + dp4data[,5])
  colnames(dp4data)[7] <- paste0(substr(x,1,4),"_Alt")
  dp4data <- cbind(dp4data, Total = dp4data[,paste0(substr(x,1,4),"_Ref")] + dp4data[, paste0(substr(x,1,4),"_Alt")])
  if(is.na(dp4datafull)){                                       # First file, use all of the SNPs
    dp4datafull <- dp4data[,c(1,6,7)]
  }else{                                                        # Otherwise we take the overlap between the existing SNPs and the new SNPs
    dp4datafull <- dp4datafull[which(dp4datafull[,1] %in% dp4data[,1]),]
    dp4data <- dp4data[which(dp4data[,1] %in% dp4datafull[,1]),]
    dp4datafull <- cbind(dp4datafull, dp4data[,c(6,7)])
  }
  cat("Number of SNps remaining:", nrow(dp4datafull), "\n")
}
write.table(dp4datafull, "allsamples.dp4",sep="\t",quote=FALSE, row.names=FALSE)
