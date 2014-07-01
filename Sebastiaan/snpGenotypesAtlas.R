# Analysis of Atlas Data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#

setwd("E:/Atlas/")

arrays <- c(#"Array1_Atlas_BioLabs_2010-02-15/BRO_Genotypes/orig/quant-norm.pm-only.brlmm-p.calls",           # I THINK THIS ONE IS NOT TO BE TRUSTED
            "Array2_Atlas_BioLabs_2010-08-20/BRO_Genotypes/orig/quant-norm.pm-only.brlmm-p.calls.txt",
            "Array3_Atlas_BioLabs_2011-10-06/BRO_Genotypes/orig/Brockmann_20111005_numcode.calls",
            "Array4_Atlas_BioLabs_2012-04-11/BRO_Genotypes/orig/quant-norm.pm-only.brlmm-p.calls",
            "Array5_Atlas_BioLabs_2012-06-25/BRO_Genotypes/orig/quant-norm.pm-only.brlmm-p.calls",
            "Array6_Atlas_BioLabs_2014-02-20/Genotypes/orig/quant-norm.pm-only.brlmm-p.calls")

calldata <- NULL; x <- 1
for(carray in arrays){
  cdata <- read.table(carray, header=TRUE, row.names=1, na.strings = "-1")            # Read in calls from a single array
  if(x == 1){ calldata <- cdata; }else{ calldata <- cbind(calldata, cdata); }
  x = x + 1
}

newdata <- calldata[, -which(colnames(calldata) == "dbSNP_RS_ID")]                    # Work on a copy of the data, remove the dbSNP column

headerIDs <- unlist(lapply(strsplit(colnames(newdata), "_"), function(x){             # Split out the correct identifiers for the arrays
  if(length(x) == 1) return(x)
  return(paste0("X", x[2]))
}))

colnames(newdata) <- headerIDs                                                        # Set the headers on the data
colnames(newdata)[38:53] <- gsub("X", "X1", colnames(newdata)[38:53])                 # Header names do not match between annotation and raw data

cat("Atlas data", nrow(newdata), "SNPs measured on", ncol(newdata), "individuals\n")

annotation <- read.table("MouseAnnotation.txt", header=TRUE)                          # Load the annotation
which(!colnames(newdata) %in% colnames(annotation))                                   # Check if we can map everything to the mouse annotation should return integer(0)

chromosomes <- c(1:19, "X", "Y", "MT")
chrAnnotationJAX <- NULL

aa <- lapply(chromosomes, function(chr){
  chrAnnotation <- read.table(paste0("JAXannotation/chr",chr,".txt"), header=FALSE, sep='\t')       # SNP / Chromosome annotation of the mouse diversity CHIP
  chrAnnotationJAX <<- rbind(chrAnnotationJAX, chrAnnotation[,c(1,9,2,3,15,12,23)])                       # Take only the annotation of interest
})
colnames(chrAnnotationJAX) <- c("JAX_ID", "Sequence", "Allele_A", "Allele_B", "dbSNP_ID", "Probe_Start", "Nucleotide_Pos")

emptySequence <- which(nchar(as.character(chrAnnotationJAX[,"Sequence"])) == 0)                     # Remove the one without a sequence
chrAnnotationJAX <- chrAnnotationJAX[-emptySequence, ]

dupEntries <- which(duplicated(as.character(chrAnnotationJAX[,"JAX_ID"])))                          # Remove the duplicate entries
chrAnnotationJAX <- chrAnnotationJAX[-dupEntries,]

positionAnnotation <- read.table("JAXblasted.txt", sep="\t")                                        # Read in the position information from Blast
colnames(positionAnnotation)[c(1,2,9,10)] <- c("JAX_ID", "Chr", "Start", "End")
positionAnnotation <- positionAnnotation[,c(1,2,9,10)]                                              # Take the columns we need from the BLAST data

duplicatedIDX <- which(duplicated(positionAnnotation[,"JAX_ID"]))                                   # Find duplicate entries
duplicatedProbes <- unique(positionAnnotation[duplicatedIDX, "JAX_ID"])                             # Find duplicate entries
multiMappingProbes <- which(positionAnnotation[, "JAX_ID"] %in% duplicatedProbes)                   # Find all probes that map more then once into the genome
positionAnnotation <- positionAnnotation[-multiMappingProbes,]                                      # Remove all probes that map more then once into the genome
cat("Removed form annotation", length(duplicatedProbes),"probes (multiple hits in the genome)\n")

orderInchrAnnotationJAX <- match(positionAnnotation[,"JAX_ID"], chrAnnotationJAX[,"JAX_ID"])        # Match up the location information from BLAST and JAX information
probeAnnotation <- cbind(positionAnnotation, chrAnnotationJAX[orderInchrAnnotationJAX,])            # Bind them together

hasAnnot <- which(rownames(newdata) %in% probeAnnotation[,"JAX_ID"])                                # What genotype data has probeAnnotation (JAX_ID)
cat(length(hasAnnot),"out of", nrow(newdata),"genotypes have no multiple hits in the genome\n")
newdata <- newdata[hasAnnot,]                                                                       # Keep only the annotated genotype data

library(biomaRt)                                                                                    # Biomart
snp.db <- useMart("snp", dataset="mmusculus_snp")                                                   # For mouse SNPs
snps <- as.character(probeAnnotation[,"dbSNP_ID"])                                                  # The list of RS_IDs we want to retrieve
snps <- snps[-which(snps=="")]                                                                      # Remove the empty ones

if(!file.exists("SNPAnnotation.txt")){
  biomartResults <- NULL
  for(x in seq(1, length(snps), 1000)){                                                             # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(snps))                                                            # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    res.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start"),                          # Use biomart to retrieve locations and reference alleles
                       filters="snp_filter", values=snps[x:xend], mart=snp.db)
    biomartResults <- rbind(biomartResults, res.biomart)
  }
  write.table(biomartResults, file="SNPAnnotation.txt", sep="\t", row.names=FALSE)
}else{                                                                                              # If the annotation file is there use it
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("SNPAnnotation.txt", sep="\t", header=TRUE)
}

orderInResults <- match(as.character(probeAnnotation[,"dbSNP_ID"]), biomartResults[,"refsnp_id"])   # Match the results from biomaRt to probeAnnotation
probeAnnotation <- cbind(probeAnnotation, biomartResults[orderInResults,])                          # Merge

SNP_Offsets <- probeAnnotation[,"Probe_Start"] - probeAnnotation[,"Nucleotide_Pos"]                 # Offset of the nucleotide interrogating the SNP

Blast_Loc <- probeAnnotation[,"Start"] + SNP_Offsets                                                # SNP positions according to Blast
probeAnnotation <- cbind(probeAnnotation, Blast_Loc)                                                # Merge

checkRSID <- apply(probeAnnotation,1,function(x){
  return(x["chrom_start"]==x["Blast_Loc"])
})

emptyData <- matrix(c("",NA,NA,NA,NA),length(which(!checkRSID)), 5, byrow=TRUE)                         # Create empty data for the probes that fail the RSID check
probeAnnotation[which(!checkRSID),c(9,12:15)] <- emptyData                                              # Zero the RSID, and locations for these ones

orderingrequested <- c("JAX_ID", "dbSNP_ID", "Chr", "Blast_Loc", "Allele_A", "Allele_B", "allele")      # Only select the columns we're interested in
probeAnnotation <- probeAnnotation[,orderingrequested]
colnames(probeAnnotation) <- c("JAX_ID", "dbSNP_ID", "Chr", "Location", "JAX_A", "JAX_B", "Allele")     # Do some renaming of columns

orderInAnnotation <- match(rownames(newdata), as.character(probeAnnotation[,"JAX_ID"]))                 # Match the annotation to the data
fulldata <- cbind(probeAnnotation[orderInAnnotation,], newdata)                                         # Bind everything together
write.table(fulldata, file="SNPAnnotated.txt", sep="\t", row.names=FALSE)                               # Write out the annotated numeric genotypes

chromosomes <- c(1:19, "X", "Y", "MT")
for(chr in chromosomes){
  cat("Number of SNPs of chr", chr, ":", length(which(fulldata[,"Chr"]==chr)), "\n")
}
