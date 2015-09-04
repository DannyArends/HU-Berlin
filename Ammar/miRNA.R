# miRNA analysis on Cows 3' UTR
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2015
# first written Aug, 2015

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

## Miranda

library(biomaRt)
mart      <- useMart("ensembl", "btaurus_gene_ensembl")
utr_seq   <- getSequence(seqType="5utr", mart=mart, type="ensembl_gene_id",id="ENSBTAG00000007695")
utr_coords <- getBM(attributes=c('ensembl_gene_id', "ensembl_transcript_id", '5_utr_start', '5_utr_end'), filters='ensembl_gene_id', values="ENSBTAG00000007695", mart=mart)
## Create the fasta file from the MIRbase reference

mdata <- read.table("miRNA.txt",sep="\t",header=TRUE)     # Data from mirbase

miRNA.bta <- mdata[which(grepl("bta", mdata[,"ID"])),]
cat("", file="cow.miRNA.fasta")

for(x in 1:nrow(miRNA.bta)){
  cat(paste0(">", miRNA.bta[x,"Mature1_ID"], "\n"), file="cow.miRNA.fasta",append=TRUE)
  cat(paste0(miRNA.bta[x,"Mature1_Seq"], "\n"), file="cow.miRNA.fasta",append=TRUE)
  if(miRNA.bta[x,"Mature2_ID"]!=""){
    cat(paste0(">", miRNA.bta[x,"Mature2_ID"], "\n"), file="cow.miRNA.fasta",append=TRUE)
    cat(paste0(miRNA.bta[x,"Mature2_Seq"], "\n"), file="cow.miRNA.fasta",append=TRUE)
  }
}

#Run prediction
execute("./miRanda-3.3a/src/miranda cow.miRNA.fasta locationsandUTR.fasta | grep '^>bta' > scanout.txt", "scanout.txt")



setwd("E:/Cow/RNA/miRNA")
SNPs <- read.table("SNPs.txt",sep="\t",header=TRUE,check.names=FALSE)       # Load the SNP data

scan.results <- read.table("scanout.txt")                                   # Load the scanresults
scan.results[,1] <-  gsub(">bta","bta",scan.results[,1])
colnames(scan.results) <- c("miRNA_id", "gene_id","score", "energy","from.miRNA","to.miRNA", "from.UTR", "to.UTR", "length", "alignment%", "alignment max%")

write.table(scan.results[which(scan.results[,"energy"] < -10),],"miRNA.scan.output",sep="\t",row.names=FALSE)   # Significant scan results (energy < -10)

significant <- scan.results[which(scan.results[,"energy"] < -10),]

genes <-  gsub(">","", readLines("locationsandUTR.fasta")[1:24 %% 2 == 1])


for(gene in genes){
  miRNAs <- scan.results[which(scan.results[,2] == gene),]
  gene_id <- strsplit(gene,"-")[[1]][2]
  prime <- strsplit(strsplit(gene,"-")[[1]][3],"")[[1]][1]
  utr.sequence <- getSequence(seqType = paste0(prime, "utr"), mart = mart, type = "ensembl_gene_id",id = gene_id)
  utr.coords <- getBM(attributes=c('ensembl_gene_id', paste0(prime,'_utr_start'), paste0(prime,'_utr_end')), filters='ensembl_gene_id', values = gene_id, mart=mart)
  
}


## PITA

#Install Pita
execute("wget http://genie.weizmann.ac.il/pubs/mir07/64bit_exe_pita_prediction.tar.gz", "64bit_exe_pita_prediction.tar.gz")
execute("tar xvfz 64bit_exe_pita_prediction.tar.gz", "Makefile")
execute("make install", "pita_prediction.pl")

#Run prediction
execute(paste0("./pita_prediction.pl -utr ", file.UTR, "-mir ", file.miRNA," -prefix ", outputname))