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

annotation <- read.table("MouseAnnotation.txt", header=TRUE)

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

which(!colnames(newdata) %in% colnames(annotation))                                   # Check if we can map everything to the mouse annotation should return integer(0)

annotation <- annotation[1:6, match(colnames(newdata), colnames(annotation))]         # Reduce the header for only the Atlas mice
rownames(annotation) <- c("Date", "Line", "Generation", "Phenotype", "Sex", "Diet")   # Correct header names
annotation <- t(annotation)                                                           # I just like annotations in columns

chromosomes <- c(1:19, "M", "X", "Y")
chrAnnotationShort <- NULL
aa <- lapply(chromosomes, function(chr){
  chrAnnotation <- read.table(paste0("JAXannotation/chr",chr,".txt"), header=FALSE, sep='\t')       # SNP / Chromosome annotation of the mouse diversity CHIP
  notperfect <- unique(as.character(chrAnnotation[which(chrAnnotation[,6] != 2), 1]))               # JAXidentifiers of not perfect hybridizing probes
  weirdProbes <- which(as.character(chrAnnotation[,1]) %in% notperfect)                             # Which ones are they
  if(length(weirdProbes) > 0){
    chrAnnotation <- chrAnnotation[-weirdProbes,]                                                   # Remove them
    cat("Chromosome", chr, "removing", length(weirdProbes), "probes\n")
  }
  chrShort <- chrAnnotation[,c(1,2,3,11,15)]                                                        # Take only the annotation of intrest
  chrAnnotationShort <<- rbind(chrAnnotationShort, chrShort)
})

noRS <- which(as.character(chrAnnotationShort[,5]) == "")                                           # Remove the ones without an RSid
chrAnnotationShort <- chrAnnotationShort[-noRS,]
dupEntries <- which(duplicated(as.character(chrAnnotationShort[,1])))                               # Remove the duplicate entries
chrAnnotationShort <- chrAnnotationShort[-dupEntries,]
colnames(chrAnnotationShort) <- c("JAX_ID", "Allele_A", "Allele_B", "Chr", "dbSNP_ID")

hasAnnot <- which(rownames(newdata) %in% chrAnnotationShort[,1])                                    # What genotype data has Annotation (dbSNP_ID)
newdata <- newdata[hasAnnot,]                                                                       # Keep only the annotated genotype data

library(biomaRt)
snp.db <- useMart("snp", dataset="mmusculus_snp")
results <- NULL
snps <- as.character(chrAnnotationShort[,"dbSNP_ID"])
for(x in seq(1, length(snps),1000)){
  xend <- min((x+1000),length(snps))
  cat("retrieving", x, "/", xend,"\n")
  res.biomart <- getBM(c("refsnp_id","allele","chr_name","chrom_start"),                                       # Use biomart to retrive the locations
                     filters="snp_filter", values=snps[x:xend], mart=snp.db)
  results <- rbind(results, res.biomart)
}

