# annotations.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014
#
# Annotating the RNA sequencing data done on our own server

library(biomaRt)

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

RPKM <- read.csv("Analysis/RPKM.txt", sep="\t", header=TRUE, check.names=FALSE)                                                          # RNA-Seq RPKM data
annotation <- read.csv("FastQ/sampledescription.txt", sep="\t", header=TRUE)                                                             # Sample annotation
RPKM <- RPKM[, -2]                                                                                                                       # Throw away the second column
colnamesShort <- unlist(lapply(strsplit(colnames(RPKM),"_"),"[",1))

colnames(RPKM) <- annotation[match(colnamesShort, annotation[,"Lib_id"]),"our_name"]                                                     # Correct column names

if(!file.exists("Analysis/BiomartAnnotation.txt")){
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")                                                # Biomart for mouse genes
  ENSMUSG <- as.character(rownames(RPKM))                                                                         # Ensemble genes we want to retrieve
  biomartResults <- NULL
  for(x in seq(0, length(ENSMUSG), 1000)){                                                                        # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(ENSMUSG))                                                                       # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")

    res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_id", "mgi_description"), 
                        filters = c("ensembl_gene_id"), 
                        values = ENSMUSG[x:xend], mart = bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
    Sys.sleep(1)
  }
  write.table(biomartResults, file="Analysis/BiomartAnnotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("Analysis/BiomartAnnotation.txt", sep="\t", header=TRUE)
}

RPKM <- RPKM[,-c(11,12)]                                                                                              # Remove the quadriceps samples

# Liver
Lcross1 <- c("F1-V-1004_L", "F1-V-1016_L", "F1-V-1020_L")
Lcross2 <- c("F1-V-1000_L", "F1-V-1008_L", "F1-V-1012_L")
LD1 <- RPKM[,Lcross1]                                                                                                 # BFMI cross BFMI860-12xB6N (D1)
LD2 <- RPKM[,Lcross2]                                                                                                 # BFMI cross B6NxBFMI860-12 (D2)
pval <- apply(cbind(LD1,LD2),1,function(x){ return(t.test(x[1:3], x[4:6])$p.value)})                                  # T-test for differences
LF1 <- cbind(LD1, LD2, apply(LD1, 1, mean), apply(LD2, 1, mean), apply(LD1, 1, mean) / apply(LD2, 1, mean), pval)
colnames(LF1)[7:10] <- c("Mean BFMI860-12xB6N L", "Mean B6NxBFMI860-12 L", "Ratio", "tTest")

LBFMIA  <- c("86026522_L", "86026502_L")
LB6NA   <- c("1006954_L", "1006956_L")
LBFMI <- RPKM[,LBFMIA]
LB6N  <- RPKM[,LB6NA]
pval <- apply(cbind(LBFMI,LB6N),1,function(x){ return(t.test(x[1:2], x[3:4])$p.value)})
LP <- cbind(LBFMI, LB6N, apply(LBFMI, 1, mean), apply(LB6N, 1, mean), apply(LBFMI, 1, mean) / apply(LB6N, 1, mean), pval)
colnames(LP)[5:8] <- c("Mean BFMI860", "Mean B6N", "Ratio", "tTest")

updatedData <- cbind(biomartResults[,c("g", "mgi_id", "symbol","chromosome_name", "start_position", "end_position", "strand", "biotype", "mgi_description")], LF1, LP)
