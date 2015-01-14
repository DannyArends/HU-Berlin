# Analysis of FV3 RNA expression data collected using Illumina microarrays
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

library(biomaRt)                                                                                                      # Biomart - Used to annotate probes
library(preprocessCore)                                                                                               # preprocessCore - For quantile normalisation

setwd("E:/Mouse/RNA/FV3")

rawdata  <- read.table("RawData/expressions.txt", sep="\t", header=TRUE, row.names=1, check.names=FALSE)              # Load the expression data
arrays      <- read.table("Annotation/arrays.txt", sep="\t", header=TRUE, row.names=1)                                # Arrays annotation
probes      <- read.table("Annotation/illumina.txt", sep="\t", header=TRUE)                                           # Probe annotation

# Create a fasta file with the probe sequences to blast against the reference genome
if(!file.exists("Analysis/probes.fasta")){
  cat("", file="Analysis/probes.fasta")
  for(x in 1:nrow(probes)){
    if(as.character(probes[x,"SEQUENCE"]) != ""){
      cat(">", probes[x,"Array_Address_Id"], "\n", sep = "", file="Analysis/probes.fasta", append=TRUE)
      cat(as.character(probes[x,"SEQUENCE"]),"\n", sep = "", file="Analysis/probes.fasta", append=TRUE)
    }
  }
}

# blastn -task blastn -query Analysis/probes.fasta -db E:/Mouse/DNA/DiversityArray/Analysis/Mus_musculus.GRCm38.74.dna.db -perc_identity 100 -outfmt 6 -evalue=0.1 -out Analysis/probelocations.txt

# Load the locations of the probes and filter them so that we only keep the unique matching probes
locations <- read.csv("Analysis/probelocations.txt", sep = "\t", header=FALSE)
colnames(locations) <- c("ProbeName", "Chr", "Ident", "Length", "U1", "U2", "U3", "Match", "Start", "Stop", "evalue", "Score")

locations <- locations[-which(locations[,"Score"] < 80),]                                                                           # Match is not good enough to be considered as duplicate ( evalue < 80 )

dupprobes <- unique(locations[which(duplicated(locations[,"ProbeName"])),"ProbeName"])                                              # Probes which have multiple matches
bestlocs <- locations[which(!duplicated(locations[,"ProbeName"])),]                                                                 # Only look up each probes once (best match)

mart      <- useMart("ensembl", "mmusculus_gene_ensembl")

# Use biomaRt to download all the genes with locations and description
if(!file.exists("Analysis/GENES.txt")){
  allgenes  <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position"), mart = mart)
  write.table(allgenes, file="Analysis/GENES.txt", sep="\t", row.names=FALSE)
}else{                                                                                           # If the annotation file is there use it
  cat("Loading biomart annotation from disk\n")
  allgenes <- read.table("Analysis/GENES.txt", sep="\t", header=TRUE)
}

# Use biomaRt to download all the exons of the genes with locations and description
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

# Match our probes to known exons within the genes
if(!file.exists("Annotation/probeannotation.txt")){
  annotationmatrix <- NULL
  for(x in 1:nrow(bestlocs)){
    inEXON <- which(as.character(allexons[,"chromosome_name"]) == as.character(bestlocs[x,"Chr"]) &                                             # In which exons ?
                    allexons[,"exon_chrom_start"] < bestlocs[x,"Start"] & allexons[,"exon_chrom_end"] > bestlocs[x,"Stop"])

    annotation <- allgenes[which(as.character(allgenes[,"ensembl_gene_id"]) %in% as.character(unique(allexons[inEXON,"ensembl_gene_id"]))), ]   # In which genes ?
    if(dim(annotation)[1] > 0){
      annotationmatrix <- rbind(annotationmatrix , cbind(ProbeName = bestlocs[x,"ProbeName"], ProbeChr = bestlocs[x,"Chr"], ProbeStart = bestlocs[x,"Start"], ProbeStop = bestlocs[x,"Stop"], annotation))
    }
    if(x %% 100 == 0){    # Print some progress report every 100 probes
      cat("Done", paste0(x,"/",nrow(bestlocs),", matched"), length(unique(annotationmatrix[,"ProbeName"])), "probes to", length(unique(annotationmatrix[,"ensembl_gene_id"])), "genes\n")
    }
  }
  cat("Matched", length(unique(annotationmatrix[,"ProbeName"])), "probes to", length(unique(annotationmatrix[,"ensembl_gene_id"])), "genes\n")  # Matched 34012 probes to 20897 genes
  write.table(annotationmatrix, file="Annotation/probeannotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading probe annotation from disk\n")
  annotationmatrix <- read.table("Annotation/probeannotation.txt", sep="\t", header=TRUE)
}

