# Analysis of Atlas Data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014

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

newdata <- calldata[, -which(colnames(calldata) == "dbSNP_RS_ID")]                    # Work on a copy of the data, remove the dbSNP column which is totaly wrong

headerIDs <- unlist(lapply(strsplit(colnames(newdata), "_"), function(x){             # Split out the correct identifiers for the arrays
  if(length(x) == 1) return(x)
  return(paste0("X", x[2]))
}))

colnames(newdata) <- headerIDs                                                        # Set the headers on the data
colnames(newdata)[38:53] <- gsub("X", "X1", colnames(newdata)[38:53])                 # Header names do not match between annotation and raw data

cat("Atlas data", nrow(newdata), "SNPs measured on", ncol(newdata), "individuals\n")

annotation <- read.table("MouseAnnotation.txt", header=TRUE)                          # Load the annotation
which(!colnames(newdata) %in% colnames(annotation))                                   # Check if we can map everything to the mouse annotation should return integer(0)

annotation <- annotation[1:6, match(colnames(newdata), colnames(annotation))]         # Reduce the header for only the Atlas mice
rownames(annotation) <- c("Date", "Line", "Generation", "Phenotype", "Sex", "Diet")   # Correct header names
annotation <- t(annotation)                                                           # I just like annotations in columns

chromosomes <- c(1:19, "X", "Y", "M")
chrAnnotationJAX <- NULL
nonUniqueCount <- 0
aa <- lapply(chromosomes, function(chr){
  chrAnnotation <- read.table(paste0("JAXannotation/chr",chr,".txt"), header=FALSE, sep='\t')       # SNP / Chromosome annotation of the mouse diversity CHIP
  notperfect <- unique(as.character(chrAnnotation[which(chrAnnotation[,6] != 2), 1]))               # JAXidentifiers of not perfect hybridizing probes
  weirdProbes <- which(as.character(chrAnnotation[,1]) %in% notperfect)                             # Which ones are they
  if(length(weirdProbes) > 0){
    chrAnnotation <- chrAnnotation[-weirdProbes,]                                                   # Remove them
    cat("Chromosome", chr, "removing", length(weirdProbes), "probes\n")
    nonUniqueCount <<- nonUniqueCount + length(weirdProbes)
  }
  chrShort <- chrAnnotation[,c(1,2,3,11,15)]                                                        # Take only the annotation of intrest
  chrAnnotationJAX <<- rbind(chrAnnotationJAX, chrShort)
})
cat("Total:", nonUniqueCount, "probes removed because of non-unique mapping\n")

dupEntries <- which(duplicated(as.character(chrAnnotationJAX[,1])))                               # Remove the duplicate entries
chrAnnotationJAX <- chrAnnotationJAX[-dupEntries,]
colnames(chrAnnotationJAX) <- c("JAX_ID", "Allele_A", "Allele_B", "Chr", "dbSNP_ID")

noRS <- which(as.character(chrAnnotationJAX[,"dbSNP_ID"]) == "")                                  # Remove the ones without an RSid
cat("Removing", length(noRS), "SNPs without an RS identifier\n")
chrAnnotationJAX <- chrAnnotationJAX[-noRS,]

hasAnnot <- which(rownames(newdata) %in% chrAnnotationJAX[,1])                                    # What genotype data has Annotation (dbSNP_ID)
newdata <- newdata[hasAnnot,]                                                                       # Keep only the annotated genotype data

library(biomaRt)                                              # Biomart
snp.db <- useMart("snp", dataset="mmusculus_snp")             # For mouse SNPs
snps <- as.character(chrAnnotationJAX[,"dbSNP_ID"])         # The RS_IDs we want to retrieve

if(!file.exists("SNPAnnotation.txt")){
  results <- NULL
  for(x in seq(1, length(snps), 1000)){                                                                    # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(snps))                                                                    # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    res.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start"),                                # Use biomart to retrieve locations and reference alleles
                       filters="snp_filter", values=snps[x:xend], mart=snp.db)
    results <- rbind(results, res.biomart)
  }
  write.table(results, file="SNPAnnotation.txt", sep="\t", row.names=FALSE)
}else{                                                                                                    # If the annotation file is there use it
  cat("Loading biomart annotation from disk\n")
  results <- read.table("SNPAnnotation.txt", sep="\t", header=TRUE)
}

orderInResults <- match(as.character(chrAnnotationJAX[,5]), results[,1])                                # Match the results from biomaRt to chrAnnotationJAX
chrAnnotationLong <- cbind(chrAnnotationJAX, results[orderInResults,])                                  # Merge

orderingrequested <- c("JAX_ID", "dbSNP_ID", "chr_name", "chrom_start", "Allele_A", "Allele_B", "allele")      # Only select the columns we're interested in
chrAnnotationLong <- chrAnnotationLong[,orderingrequested]
colnames(chrAnnotationLong) <- c("JAX_ID", "dbSNP_ID", "Chr", "Location", "JAX_A", "JAX_B", "Allele")     # Do some renaming of columns

orderInAnnotation <- match(rownames(newdata), as.character(chrAnnotationLong[,"JAX_ID"]))                 # Match the annotation to the data
fulldata <- cbind(chrAnnotationLong[orderInAnnotation,], newdata)                                         # Bind everything together
write.table(fulldata, file="SNPAnnotated.txt", sep="\t", row.names=FALSE)                                 # Write out the annotated numeric genotypes

chromosomes <- c(1:19, "X", "Y", "M")
for(chr in chromosomes){
  cat("Number of SNPs of chr",chr,":",length(which(fulldata[,"Chr"]==chr)),"\n")
}
