# miRNA analysis on Cows 3' UTR
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2015
# first written Aug, 2015

## Miranda

setwd("E:/Cow/RNA/miRNA")
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

execute("./miRanda-3.3a/src/miranda cow.miRNA.fasta reference.txt | grep '^>bta' > scanout.txt", "scanout.txt")

scan.results <- read.table("scanout.txt")
scan.results[,1] <-  gsub(">bta","bta",scan.results[,1])
colnames(scan.results) <- c("miRNA_id", "gene_id","score", "energy","from.miRNA","to.miRNA", "from.UTR", "to.UTR", "length", "alignment%", "alignment max%")

write.table(scan.results[which(scan.results[,"energy"] < -10),],"miRNA.scan.output",sep="\t",row.names=FALSE)


## PITA

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

execute("wget http://genie.weizmann.ac.il/pubs/mir07/64bit_exe_pita_prediction.tar.gz", "64bit_exe_pita_prediction.tar.gz")
execute("tar xvfz 64bit_exe_pita_prediction.tar.gz", "Makefile")
execute("make install", "pita_prediction.pl")

execute(paste0("./pita_prediction.pl -utr ", file.UTR, "-mir ", file.miRNA," -prefix ", outputname))