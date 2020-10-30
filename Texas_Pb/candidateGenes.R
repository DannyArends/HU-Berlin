setwd("D:/Edrive/Mouse/Texas_Pb")

library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

regions <- c("1:124548738:184926264", 
             "1:135464050:157672910", 
             "3:85146443:147514132", 
             "4:58831312:88222686", 
             "6:76940963:114257130", 
             "7:62305639:81923457", 
             "7:62305639:75768444", 
             "7:19482853:119720011", 
             "8:31799796:104944836")
for(r in regions){
  res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                          filters = c("chromosomal_region", "biotype"), values = list(r, "protein_coding"), mart = bio.mart)
  write.table(res.biomart, file=paste0("genes_", gsub(":", "-",r), ".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}