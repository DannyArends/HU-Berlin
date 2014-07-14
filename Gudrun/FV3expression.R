# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014

setwd("E:/Mouse/RNA/FV3")

library(biomaRt)                                                                                                      # Biomart - Used to annotate probes
library(preprocessCore)                                                                                               # preprocessCore - For quantile normalisation

expression  <- read.table("FV3_expressionData_unnormalized.txt", sep="\t", header=TRUE, row.names=1)                  # Load the expression data

probes      <- read.table("probeannotation.txt", sep="\t", header=TRUE)                                               # Probe annotation, keep only the Array_Address_Id and the Sequence
probes      <- probes[, c("Array_Address_Id", "SEQUENCE")]

arrays      <- read.table("arrayannotation.txt", sep="\t", header=TRUE, row.names=1)                                  # Arrays annotation
rownames(arrays) <- paste0("X",rownames(arrays))                                                                      # Fix the rownames, such that they match the column names of expression data

# Pre-processing steps to create a fast sequence file we can use to BLAST against the mouse genome
if(!file.exists("FV3sequences.fasta")){
  cat("", file="FV3sequences.fasta")                                                                                  # Create a fasta file to blast the probe sequences to the mouse genome
  aa <- apply(probes, 1, function(x){
    cat(paste0(">", as.numeric(x[1]), "\n", x[2], "\n"), file="FV3sequences.fasta", append=TRUE)
  })
}

# After this use blastn to query the database for the location of the JAX probe sequences
# blastn -task blastn -query FV3sequences.fasta -db E:/Mouse/DNA/DiversityArray/Analysis/Mus_musculus.GRCm38.74.dna.db -perc_identity 100 -outfmt "6 qseqid sseqid sstart send sstrand" -evalue 0.1 -num_alignments 5 -out IlluminaProbesBlast.txt

# Pre-processing data, Log transform, and normalize match the array annotations

expdata <- apply(expression, 2, as.numeric)                                                                           # Force a numeric matrix
colnames(expdata) <- colnames(expression); rownames(expdata) <- rownames(expression)                                  # Add the column and row names

expdata <- log(expdata)                                                                                               # Log transform
boxplot(expdata, las=2)

expdata <- expdata[, -which(colnames(expdata) == "X1740174012_D")]                                                    # Misbehaving array X1740174012_D, expression = 10 x lower the expected
boxplot(expdata, las=2)

normdata <- normalize.quantiles(expdata)                                                                              # Basic Quantile normalisation
boxplot(normdata, las=2)
colnames(normdata) <- colnames(expdata); rownames(normdata) <- rownames(expdata)

hasAnnotation <- match(colnames(normdata), rownames(arrays))                                                          # Annotation matching to normdata
arrays <- arrays[hasAnnotation, ]

# Biomart analysis of the BLAST results

blastResults <- read.table("IlluminaProbesBlast.txt", sep="\t", colClasses="character")                               # Load out custom blast results
multiMatch <- blastResults[duplicated(blastResults[,1]),1]                                                            # Probes with more then 1 unique hit in the genome
blastResults <- blastResults[-which(blastResults[,1] %in% multiMatch), ]                                              # Remove them from the data

blastResults <- cbind(blastResults, NA)                                                                               # Additional column for the gene annotation
colnames(blastResults) <- c("decoder_ID", "Chr", "Start", "Stop", "Strand", "Annotation")                             # Name the columns

blastResults[which(blastResults[,"Strand"] == "plus"),"Strand"]  <- 1
blastResults[which(blastResults[,"Strand"] == "minus"),"Strand"] <- -1
chr.regions <- paste0(blastResults[,"Chr"], ":",blastResults[,"Start"],":",blastResults[,"Stop"])                     # Create the chromosome regions we want to query

if(!file.exists("BiomartAnnotation.txt")){
  bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
  biomartResults <- NULL
  for(x in seq(0, length(chr.regions), 1000)){                                                                        # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(chr.regions))                                                                       # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")

    res.biomart <- getBM(attributes = c("mgi_id","mgi_symbol", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(chr.regions[x:xend],"protein_coding"), mart = bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
    Sys.sleep(1)
  }
  write.table(biomartResults, file="BiomartAnnotation.txt", sep="\t", row.names=FALSE)
}else{                                                                                                                # Or load them from Disk if we already have them
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("BiomartAnnotation.txt",sep="\t", header=TRUE)
}

biomartResults <- biomartResults[-which(duplicated(biomartResults[,"mgi_id"])),]                                      # Remove the duplicate annotations

# Merge the biomart analysis to the BLAST results
if(!file.exists("IlluminaProbesBlastAnnotated.txt")){
  for(x in 1:nrow(blastResults)){
    mgis <- unique(biomartResults[which(biomartResults[,"chromosome_name"] == blastResults[x,"Chr"] & 
                  biomartResults[,"start_position"] <= as.numeric(blastResults[x,"Start"]) & 
                  biomartResults[,"end_position"] >= as.numeric(blastResults[x,"Stop"]) &
                  biomartResults[,"strand"] == blastResults[x,"Strand"]
                  ),"mgi_id"])
    if(length(mgis) > 0){
      annot <- NULL
      for(hit in mgis){
        id <- which(biomartResults[,"mgi_id"] == hit)
        annot <- c(annot, paste(hit, biomartResults[id,"mgi_symbol"], biomartResults[id,"start_position"], biomartResults[id,"end_position"],sep=";"))
      }
      blastResults[x,"Annotation"] = paste0(annot,collapse="///")
    }
    if(x %% 1000 == 0) cat("Done",x,"/",nrow(blastResults),"\n")
  }
  write.table(blastResults, file="IlluminaProbesBlastAnnotated.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading annotated BLAST results from disk\n")
  blastResults <- read.table("IlluminaProbesBlastAnnotated.txt",sep="\t", header=TRUE, colClasses=c("character"))
}

# Not expressed genes
expressionlevels <- cbind(apply(normdata, 1, median), apply(normdata, 1, sd))                             # Calculate the expression levels

notExpressed <- which(expressionlevels[,1] < 4.65 & expressionlevels[,2] < 0.1)                            # Mean expression lower then 4.6 and SD < 0.1
image(normdata[notExpressed,])

cat("Found", length(notExpressed), "Genes/Probes\n")

hasAnnotation <- match(rownames(normdata[notExpressed,]), blastResults[,"decoder_ID"])
hasAnnotation <- na.omit(hasAnnotation)

notExpressed <- na.omit(blastResults[hasAnnotation,])

notExpressedAnnot <- lapply(unlist(strsplit(notExpressed[,"Annotation"],"///")), strsplit, ";")
notExpressedMGI <- unlist(lapply(notExpressedAnnot,function(x){return(x[[1]][1]);}))

cat("Not expressed MGIs", length(notExpressedMGI),"before using RNASeq data\n")

setwd("E:/Mouse/RNA/Sequencing")                                                                          # Adding a second data source RNA Seq expression data

RNASeqData <- read.table("BFMI_RPKM_ANN.txt", sep="\t", header=TRUE)

RNASeqData <- RNASeqData[which(RNASeqData[,"mgi_id"] %in% notExpressedMGI),]
MeanRNASeq <- RNASeqData[,grep("Mean", colnames(RNASeqData))]

RNASeqExpressedMGI <- RNASeqData[which(!apply(MeanRNASeq,1,mean) < 15), "mgi_id"]                         # Expressed in the RNASeq dataset so not remove them

notExpressedMGI <- notExpressedMGI[-which(notExpressedMGI %in% RNASeqExpressedMGI)]
cat("Not expressed MGIs", length(notExpressedMGI),"after using RNASeq data\n")

setwd("E:/Mouse/RNA/FV3")                                                                                 # Return to the FV3 Dataset

cat(notExpressedMGI, sep="\n",file="notExpressedMGIs.txt")

# Use http://www.informatics.jax.org/batch to query the resulting to REFSEQ sequences, then continue in the folder Array Design

# Differential expressions of genes: tissue, line, diet and line x diet

celltypes <- c("fat", "liver", "brain")
condition <- c("FF", "NF")

pvalues <- matrix(NA, nrow(normdata), 4)
for(x in 1:nrow(normdata)){
  pvalues[x, 1:4] <- anova(lm(normdata[x,] ~ arrays[,"tissue"] + arrays[,"line"] * arrays[,"diet"]))[[5]][1:4]
  if(x %% 1000 == 0) cat("Done", x, "/", nrow(normdata), "\n")
}
