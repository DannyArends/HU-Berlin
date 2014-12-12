library(biomaRt)
setwd("E:/Mouse/RNA/ArrayDesign/Atlas data")

mart      <- useMart("ensembl", "mmusculus_gene_ensembl")
allgenes  <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), mart = mart)

if(!file.exists("Analysis/EXONS.txt")){
  allexons <- NULL
  for(x in seq(1, length(allgenes[,"ensembl_gene_id"]), 1000)){                                  # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(allgenes[,"ensembl_gene_id"]))                                 # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    exons  <- getBM(attributes = c("ensembl_gene_id", "ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end"), filter="ensembl_gene_id", values=allgenes[,"ensembl_gene_id"][x:xend], mart = mart)
    allexons <- rbind(allexons, exons)
  }
  write.table(allexons, file="Analysis/EXONS.txt", sep="\t", row.names=FALSE)
}else{                                                                                           # If the annotation file is there use it
  cat("Loading biomart annotation from disk\n")
  allexons <- read.table("Analysis/EXONS.txt", sep="\t", header=TRUE)
}




