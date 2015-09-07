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
mdata       <- read.table("miRNA.txt",sep="\t",header=TRUE)     # Data from mirbase
miRNA.bta   <- mdata[which(grepl("bta", mdata[,"ID"])),]        # Take the cow known miRNA

cat("", file="cow.miRNA.fasta")                                 # Empty the output fasta, and fill it
for(x in 1:nrow(miRNA.bta)){
  cat(paste0(">", miRNA.bta[x,"Mature1_ID"], "\n"), file="cow.miRNA.fasta",append=TRUE)
  cat(paste0(miRNA.bta[x,"Mature1_Seq"], "\n"), file="cow.miRNA.fasta",append=TRUE)
  if(miRNA.bta[x,"Mature2_ID"]!=""){
    cat(paste0(">", miRNA.bta[x,"Mature2_ID"], "\n"), file="cow.miRNA.fasta",append=TRUE)
    cat(paste0(miRNA.bta[x,"Mature2_Seq"], "\n"), file="cow.miRNA.fasta",append=TRUE)
  }
}

#Run prediction (on server, since miranda is a linux tool)
execute("./miRanda-3.3a/src/miranda cow.miRNA.fasta locationsandUTR.fasta | grep '^>bta' > scanout.txt", "scanout.txt")

# Transfer the scanout back to the local machine
#library(biomaRt)
#mart      <- useMart("ensembl", "btaurus_gene_ensembl")
setwd("E:/Cow/RNA/miRNA")

SNPfasta <- read.table("Something.txt")
seqs <- SNPfasta[1:nrow(SNPfasta) %% 2 == 0,]
names(seqs) <- SNPfasta[1:nrow(SNPfasta) %% 2 == 1,]
names(seqs) <-  gsub(">","", names(seqs))    # get the gene names from the FASTA (we do not use the location)


scan.results <- read.table("scanout.txt",colClasses="character")                                     # Load the scanresults
scan.results[,1] <-  gsub(">bta","bta",scan.results[,1])
colnames(scan.results) <- c("miRNA_id", "gene_id","score", "energy","from.miRNA","to.miRNA", "from.UTR", "to.UTR", "length", "alignment%", "alignment max%")

significant <- scan.results[which(as.numeric(scan.results[,"energy"]) < -10),]            # Filter the results for a maximum binding energy of -10
write.table(significant, "miRNA.scan.output", sep="\t", row.names = FALSE)    # Significant scan results (energy < -10)

results <- NULL

for(name in names(seqs)){
  onUTR <- significant[which(significant[,"gene_id"] == name), ]
  mletters <- strsplit(as.character(seqs[name]),"")[[1]]
  positions <- which(!is.na(as.numeric(mletters)))
  cnt <- 1
  for(pos in positions){
    for(utr in 1:nrow(onUTR)){
      if(as.numeric(onUTR[utr,"from.UTR"]) <= pos && as.numeric(onUTR[utr,"to.UTR"]) >= pos){
        #cat("Match", name, "",cnt, onUTR[utr,]),"\n")
        results <- rbind(results, c(name, cnt, pos, as.character(onUTR[utr,])))
      }
    }
    cnt <- cnt + 1
  }
}
write.table(results,"results.txt",sep="\t")


## PITA

#Install Pita
execute("wget http://genie.weizmann.ac.il/pubs/mir07/64bit_exe_pita_prediction.tar.gz", "64bit_exe_pita_prediction.tar.gz")
execute("tar xvfz 64bit_exe_pita_prediction.tar.gz", "Makefile")
execute("make install", "pita_prediction.pl")

#Run prediction
execute(paste0("./pita_prediction.pl -utr ", file.UTR, "-mir ", file.miRNA," -prefix ", outputname))