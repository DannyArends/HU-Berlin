library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

ensemblIDS <- c("ENSMUSG00000000031")

res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("ensembl_gene_id"), 
                        values = ensemblIDS, mart = bio.mart)
                        
mgiSymblos <- c("H19")

res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("mgi_symbol"), 
                        values = mgiSymblos, mart = bio.mart)