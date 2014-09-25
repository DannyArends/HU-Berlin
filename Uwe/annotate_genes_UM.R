# RHB/UM 2014-09
#Gesprächsnotizen, Anmerkungen:
#Die listAttributes/Filter() Aufrufe sind hier auskommentiert, da diese nur zur Informationsbeschaffung 
#für den Programmierer dienen. Das spart diesen ellenlangen Output.
#Zuerst die Funktionen im unteren Teil des Skriptes einzulesen.
#Bei der for-Schleife über die SNP-Liste wird eine Funktion aufgerufen, welche die Gen-Listen in separate Dateien schreibt. 
#Dabei gibt Windows eine Warnung aus, weil wiederholt in die gleiche Datei geschrieben wird (Headerzeilen, Tabelle). 
#Diese Warnung kann man mit 1x "Cancel" klicken abstellen.
#Die resultierenden "Bt_candidate_SNPs_ARS-BFGL-NGS-XYZ.txt" Dateien am besten nach Excel importieren - ein Tabellenblatt 
#per SNP (copy/paste reicht), damit die langen Einträge der Genfunktionen/GO-Annotation, optisch abgeschnitten werden 
#und die Tabelle kompakt und übersichtlich dargestellt wird. 

# Script creates PED file with up-to-date genomic positions (as of the assembly# used in Ensembl) of user-provided SNPs (rsIDs required!)
# Script fetches genes+their annotation in a specified genomic region ("window") around the SNP positions

# ---- HELPER FUNCTIONS ----
concat_unique=function(x){ as.character(paste(sort(unique(x)),collapse=",")) }

# pos: vector with gene start positions
# refpos: reference position
# return: list with difference between position and reference position
calc_distance <- function(pos, refpos){ return(abs(as.numeric(pos)-as.numeric(refpos))); }

# data: data frame or table to be written to file
write_table_with_header = function(data, infile, outfile, info){
  write(file=outfile, "# AUTHOR: UM - HU Berlin")
  write(file=outfile, paste("#   DATE:", date()), append=TRUE)
  write(file=outfile, paste("# SCRIPT:", script), append=TRUE)
  write(file=outfile, paste("# INFILE:", infile ), append=TRUE)
  write(file=outfile, paste("# OUTFILE:", outfile), append=TRUE)
  write.table(file=outfile, data, quote=FALSE, row.names=FALSE, sep="\t", append=TRUE)
}

# ---- SCRIPT STARTS HERE ----
script <- "annotate_genes_in_SNP_region.R"                                                        # Name of the script
setwd("D:/Halde/RBB_u_ExkursionRBB/DSN_DONOREN_Gassan/Ergebnisse/2014_09/R2014_09")               # Working directory

#source("http://bioconductor.org/biocLite.R")  # execute only once for installation!
#biocLite("biomaRt")                           # execute only once for installation! 
library(biomaRt)                                                                                  # installed via Bioconductor!
df.ens.version <- listMarts()                                                                     # which database version are we dealing with?
snp.version <- as.character(df.ens.version[which(df.ens.version[,"biomart"]=="snp"), "version"])  # Determine Ensembl Variation Database Version

# User list of SNPs with rsIDs local submitter IDs (perhaps read in as table)
user.dat <- data.frame(refsnp_id = c("rs110573103",
                                     "rs109421300",
                                     "rs41657725",
                                     "rs43479230",
                                     "rs468982120",
                                     "rs466267582"),
                       snp.loc.ID = c("ARS-BFGL-NGS-31468",
                                      "ARS-BFGL-NGS-4939",
                                      "BTA-90692-no-rs",
                                      "BTB-00272812",
                                      "Hapmap39651-BTA-42671",
                                      "HapMap47255-BTA-34035"), stringsAsFactors=FALSE)

str(user.dat)

# ------------- ANNOTATE GENOMIC POSITIONS ------------------

snpmart <- useMart("snp", dataset = "btaurus_snp")                                    # Create biomaRt adapter for Ensembl SNP (variation) database

# Get SNP infos
BMattributes  <- c("refsnp_id", "allele", "chr_name", "chrom_start", "chrom_strand", "ensembl_type", "consequence_type_tv","distance_to_transcript")
snp.rsID.info <- getBM(BMattributes, filters = c("snp_filter"), values = user.dat[,"refsnp_id"], mart = snpmart)
user.dat.pos  <- merge(user.dat, snp.rsID.info, by="refsnp_id", sort=FALSE)           # data.frame of with current genomic positions annotated list of user SNPs

# Danny Arends: This is a check to see if there are missing SNPs
notFound <- (!user.dat[,"refsnp_id"] %in% user.dat.pos[,"refsnp_id"])
if(any(notFound)){ 
  cat("SEVERE WARNING: Some SNPs were not found in the database:", user.dat[notFound,"refsnp_id"],"\n");
  user.dat.pos <- rbind(user.dat.pos,                                                 # Add the missing SNPs manually, this is not a lasting solution
  c("rs468982120", "Hapmap39651-BTA-42671", "T/?", 22, 42671, 1,NA),                  # http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=rs468982120
  c("rs466267582", "HapMap47255-BTA-34035", "C/?", 17, 34035, 1,NA))                  # http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?searchType=adhoc_search&type=rs&rs=rs466267582
  user.dat.pos <- as.data.frame(user.dat.pos)
}
# When we get the SEVER WARNING, we have to manually add the locations of the SNPs (see above)

write_table_with_header(user.dat.pos, NA, "Bt_candidate_SNPs.txt", NA)                # Optionally write the returned location data to outfile

# ------------- ANNOTATE SNP ADJACENT GENES ------------------
k <- 10^6                                                                             # Define genomic region of length 2*k around SNP positions 1 Mb, i.d. 2Mb region scanned for genes potentially linked to the "center" SNP
bt.mart <- useMart("ensembl", dataset="btaurus_gene_ensembl")                         # Create specific biomaRt adapter for Ensembl cattle gene database

for(i in 1:nrow(user.dat.pos)){  # loop over user SNPS
    snp.ID <- user.dat.pos[i,"snp.loc.ID"]
    user.filter <- paste(c(
    user.dat.pos[i,"chr_name"],
    as.numeric(user.dat.pos[i,"chrom_start"])-k,
    as.numeric(user.dat.pos[i,"chrom_start"])+k,
    user.dat.pos[i,"chrom_strand"]), collapse=":")

    # Get gene information
    bt.snp.genes <- getBM( c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", "end_position", "gene_biotype", "description"), filters = "chromosomal_region", values = user.filter, mart = bt.mart)
    #cat("Found:", nrow(bt.snp.genes),"genes in range of", snp.ID,"\n")
    
    distance <- unlist(lapply(bt.snp.genes[,c("start_position")], calc_distance, user.dat.pos[1,"chrom_start"]))  # Compute distance between (center)SNP and gene starts
    bt.snp.genes.dist <- data.frame(bt.snp.genes, DIST_TO_SNP = distance, stringsAsFactors=FALSE)

    # Get GOSlim (Gene Ontology) annotation
    bt.snp.genes.goa <- getBM( c("ensembl_gene_id", "goslim_goa_description"), filters = "ensembl_gene_id", values = bt.snp.genes[,"ensembl_gene_id"], mart = bt.mart)

    # Remove gene_id redundancy by aggregating all GO annotations per gene into one string
    bt.snp.genes.goa.agg <- aggregate( x=list(GOA=bt.snp.genes.goa$goslim_goa_description), by=list(ensembl_gene_id=bt.snp.genes.goa$ensembl_gene_id), concat_unique)
    bt.snp.genes.annot   <- merge(bt.snp.genes.dist, bt.snp.genes.goa.agg, by="ensembl_gene_id", sort=F)

    # EXAMPLE: query infos about genes in the (putative!) SNP linkage region
    #bt.snp.genes.annot[grep("signal transduction",bt.snp.genes.annot[,"GOA"]), c("external_gene_id"]
    
    # Optionally write to ONE file per SNP
    cat("Found:", nrow(bt.snp.genes.annot),"annotated genes in range of", snp.ID,"\n")
    write_table_with_header(bt.snp.genes.annot, NA, paste("Bt_candidate_SNPs_", snp.ID, ".txt", sep=""), NA)
}
cat("Finished analysis of", nrow(user.dat.pos), "SNPs\n")
