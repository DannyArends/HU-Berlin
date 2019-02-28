source("https://bioconductor.org/biocLite.R")
biocLite("biomaRt")

library(biomaRt)

setwd("D:/Collegues/Monika")
dat <- read.csv("Eca-GOA.txt", sep="\t")

#57                                                 protein_id                              Protein (Genbank) ID [e.g. AAA02487]
#58                                                     refseq_mrna                                   RefSeq mRNA [e.g. NM_001195597]
#59                                           refseq_mrna_predicted                         RefSeq mRNA predicted [e.g. XM_001125684]
#60                                                    refseq_ncrna                                     RefSeq ncRNA [e.g. NR_002834]
#61                                                  refseq_peptide                             RefSeq Protein ID [e.g. NP_001005353]
#62                                        refseq_peptide_predicted                   RefSeq Predicted Protein ID [e.g. XP_001720922]

mart      <- useMart("ensembl", "ecaballus_gene_ensembl")


dat[which(grepl("XM", dat[,"Name"])),"Name"]

values <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(dat[, "Name"]), "|",fixed=TRUE),"[", 4)), ".", fixed=TRUE), "[",1))

r1 <- getBM(attributes = c("ensembl_gene_id", "refseq_mrna_predicted", "chromosome_name", "start_position", "end_position", "external_gene_name", "wikigene_description", "goslim_goa_accession", "goslim_goa_description"), filter="refseq_mrna_predicted", values = values, mart = mart)
r2 <- getBM(attributes = c("ensembl_gene_id", "refseq_mrna", "chromosome_name", "start_position", "end_position", "external_gene_name", "wikigene_description", "goslim_goa_accession", "goslim_goa_description"), filter="refseq_mrna", values = values, mart = mart)
r3 <- getBM(attributes = c("ensembl_gene_id", "refseq_ncrna", "chromosome_name", "start_position", "end_position", "external_gene_name", "wikigene_description", "goslim_goa_accession", "goslim_goa_description"), filter="refseq_ncrna", values = values, mart = mart)

colnames(r1)[2] <- "value"
colnames(r2)[2] <- "value"
colnames(r3)[2] <- "value"

results <- rbind(r1,r2,r3)

idx <- which(results[, "goslim_goa_description"] %in% c("molecular_function", "biological_process", "cellular_component"))
results <- results[-idx,]

dat <- cbind(dat, GOA...IDS = NA)
dat <- cbind(dat, external_gene_name = NA)

for(x in 1:nrow(results)){
  idx <- which(values == results[x, "value"])
  if(!is.na(dat[idx, "GOA...Function"])){
    dat[idx, "GOA...Function"] <- paste(dat[idx, "GOA...Function"], results[x, "goslim_goa_description"], sep=", ")
  }else{
    dat[idx, "GOA...Function"] <- results[x, "goslim_goa_description"]
  }
  if(!is.na(dat[idx, "GOA...IDS"])){
    dat[idx, "GOA...IDS"] <- paste(dat[idx, "GOA...IDS"], results[x, "goslim_goa_accession"], sep=", ")
  }else{
    dat[idx, "GOA...IDS"] <- results[x, "goslim_goa_accession"]
  }
  if(is.na(dat[idx, "external_gene_name"])){
    dat[idx, "external_gene_name"] <- results[x, "external_gene_name"]
  }
}

write.table(dat, file="Eca-GOA-DA.txt", sep="\t", row.names=FALSE,na="")