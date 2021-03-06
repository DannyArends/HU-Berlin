# Analysis of the micro array data from Atlas 2014
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

library(preprocessCore)
library(biomaRt)

# Loading data
setwd("D:/Edrive/Mouse/RNA/ArrayDesign/Atlas data")
arrays <- read.table("Annotation/arrays.txt", header=TRUE, sep="\t", colClasses="character")

if(!file.exists("Analysis/arraydata.txt")){
  cat("Loading raw data from disk\n")
  rawdata <- NULL
  for(filename in arrays[,"Filename"]){
    marray <- read.csv(paste0("Analyseergebnisse/",filename), sep = "\t", skip = 9, header = TRUE)
    rawdata <- cbind(rawdata, marray$gProcessedSignal)
  }
  rawdata <- cbind(marray[,c("ProbeName", "Sequence")], rawdata)
  colnames(rawdata) <- c("ProbeName", "Sequence", arrays[,"AtlasID"])

  #Some probes are on the array multiple times, so let's de-duplicate them (we take the average expression over the similar probes)
  uniqueprobes <- unique(rawdata[,"ProbeName"])
  arraydata <- matrix(NA,length(uniqueprobes), (ncol(rawdata)-1), dimnames = list(uniqueprobes, colnames(rawdata)[-1]))
  cat("De-duplicating probes\n")
  for(probe in uniqueprobes){
    prows <- which(rawdata[,"ProbeName"] == probe)
    probeseq <- rawdata[prows,"Sequence"][1]
    arraydata[probe, ] <- c(as.character(probeseq), apply(rawdata[which(rawdata[,"ProbeName"] == probe), arrays[,"AtlasID"]], 2, mean))
  }
  arraydata <- data.frame(arraydata)
  write.table(arraydata, file="Analysis/arraydata.txt", sep="\t", row.names=TRUE)
}else{
  arraydata <- read.table("Analysis/arraydata.txt", sep="\t", header=TRUE, row.names=1)
}

# Preprocessing
arraydata[,arrays[,"AtlasID"]] <- log2(arraydata[,arrays[,"AtlasID"]])                              # Log2 transformation
boxplot(arraydata[,arrays[,"AtlasID"]], las=2)
arraydata[,arrays[,"AtlasID"]] <- normalize.quantiles(as.matrix(arraydata[,arrays[,"AtlasID"]]))    # Quantile normalisation
boxplot(arraydata[,arrays[,"AtlasID"]], las=2)

# QC of our probes
boxplot(arraydata[,arrays[,"AtlasID"]], las=2)
heatmap(cor(arraydata[,arrays[,"AtlasID"]], method="spearman"))

# Create a fasta file with the probe sequences to blast against the reference genome
if(!file.exists("Analysis/probes.fasta")){
  cat("", file="Analysis/probes.fasta")
  for(x in 1:nrow(arraydata)){
    if(as.character(arraydata[x,"Sequence"]) != ""){
      cat(">", rownames(arraydata)[x], "\n", sep = "", file="Analysis/probes.fasta", append=TRUE)
      cat(as.character(arraydata[x,"Sequence"]),"\n", sep = "", file="Analysis/probes.fasta", append=TRUE)
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
    inEXON <- which(as.character(allexons[,"chromosome_name"]) == as.character(bestlocs[x,"Chr"]) &                                            # In which exons ?
                    allexons[,"exon_chrom_start"] < bestlocs[x,"Start"] & allexons[,"exon_chrom_end"] > bestlocs[x,"Stop"])

    annotation <- allgenes[which(as.character(allgenes[,"ensembl_gene_id"]) %in% as.character(unique(allexons[inEXON,"ensembl_gene_id"]))), ]   # In which genes ?
    if(dim(annotation)[1] > 0){
      annotationmatrix <- rbind(annotationmatrix , cbind(ProbeName = bestlocs[x,"ProbeName"], ProbeChr = bestlocs[x,"Chr"], ProbeStart = bestlocs[x,"Start"], ProbeStop = bestlocs[x,"Stop"], annotation))
    }
    if(x %% 100 == 0){    # Print some progress report every 100 probes
      cat("Done", paste0(x,"/",nrow(bestlocs),", matched"), length(unique(annotationmatrix[,"ProbeName"])), "probes to", length(unique(annotationmatrix[,"ensembl_gene_id"])), "genes\n")
    }
  }
  cat("Matched", length(unique(annotationmatrix[,"ProbeName"])), "probes to", length(unique(annotationmatrix[,"ensembl_gene_id"])), "genes\n")
  write.table(annotationmatrix, file="Annotation/probeannotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading probe annotation from disk\n")
  annotationmatrix <- read.table("Annotation/probeannotation.txt", sep="\t", header=TRUE)
}

# Per probe differential expression analysis
pvalues <- t(apply(arraydata[,arrays[,"AtlasID"]], 1, function(x){
  res <- anova(lm(as.numeric(x) ~ arrays[,"Tissue"] + arrays[,"Strain"]))
  return(res[[5]][1:2])
}))

colnames(pvalues) <- c("Tissue_P", "Strain_P")
arraydata         <- cbind(arraydata, pvalues)

# Filter only probes that match genes in the genome
arraydata <- arraydata[which(rownames(arraydata) %in% annotationmatrix[,"ProbeName"]),]

# Filter for differentially expressed probes
DEthreshold <- 0.1 / nrow(arraydata)
DEstrain    <- which(as.numeric(arraydata[,"Strain_P"]) < DEthreshold)
DEtissue    <- which(as.numeric(arraydata[,"Tissue_P"]) < DEthreshold)

cat("Found", length(DEtissue), "probes differentially expressed between tissues\n")
cat("Found", length(DEstrain), "probes differentially expressed between strains\n")

# The different groups in the data
HT      <- arrays[which(arrays[,"Tissue"] == "HT"),"AtlasID"]
HT_BFMI <- arrays[which(arrays[,"Tissue"] == "HT" & arrays[,"Strain"] == "BFMI"),"AtlasID"]
HT_B6N  <- arrays[which(arrays[,"Tissue"] == "HT" & arrays[,"Strain"] == "B6N") ,"AtlasID"]
HT_F1   <- arrays[which(arrays[,"Tissue"] == "HT" & arrays[,"Strain"] == "F1")  ,"AtlasID"]
GF      <- arrays[which(arrays[,"Tissue"] == "GF"),"AtlasID"]
GF_BFMI <- arrays[which(arrays[,"Tissue"] == "GF" & arrays[,"Strain"] == "BFMI"),"AtlasID"]
GF_B6N  <- arrays[which(arrays[,"Tissue"] == "GF" & arrays[,"Strain"] == "B6N") ,"AtlasID"]
GF_F1   <- arrays[which(arrays[,"Tissue"] == "GF" & arrays[,"Strain"] == "F1")  ,"AtlasID"]

# Add the mean values per group
arraydata <- cbind(arraydata, HT_Mean         = round(apply(arraydata[,HT],1,mean), 2))
arraydata <- cbind(arraydata, HT_SD           = round(apply(arraydata[,HT],1,sd), 2))
arraydata <- cbind(arraydata, GF_Mean         = round(apply(arraydata[,GF],1,mean), 2))
arraydata <- cbind(arraydata, GF_SD           = round(apply(arraydata[,GF],1,sd), 2))
arraydata <- cbind(arraydata, ratio_HT_vs_GF  = round(apply(arraydata[,HT],1,mean) / apply(arraydata[,GF],1,mean), 2))

arraydata <- cbind(arraydata, HT_BFMI_Mean    = round(apply(arraydata[,HT_BFMI],1,mean), 2))
arraydata <- cbind(arraydata, HT_BFMI_SD      = round(apply(arraydata[,HT_BFMI],1,sd), 2))
arraydata <- cbind(arraydata, GF_BFMI_Mean    = round(apply(arraydata[,GF_BFMI],1,mean), 2))
arraydata <- cbind(arraydata, GF_BFMI_SD      = round(apply(arraydata[,GF_BFMI],1,sd), 2))
arraydata <- cbind(arraydata, ratio_HT_BFMI_vs_GF_BFMI  = round(apply(arraydata[,HT_BFMI],1,mean) / apply(arraydata[,GF_BFMI],1,mean), 2))

arraydata <- cbind(arraydata, HT_B6N_Mean     = round(apply(arraydata[,HT_B6N],1,mean), 2))
arraydata <- cbind(arraydata, HT_B6N_SD       = round(apply(arraydata[,HT_B6N],1,sd), 2))
arraydata <- cbind(arraydata, GF_B6N_Mean     = round(apply(arraydata[,GF_B6N],1,mean), 2))
arraydata <- cbind(arraydata, GF_B6N_SD       = round(apply(arraydata[,GF_B6N],1,sd), 2))
arraydata <- cbind(arraydata, ratio_HT_B6N_vs_GF_B6N  = round(apply(arraydata[,HT_B6N],1,mean) / apply(arraydata[,GF_B6N],1,mean), 2))

arraydata <- cbind(arraydata, HT_F1_Mean      = round(apply(arraydata[,HT_F1],1,mean), 2))
arraydata <- cbind(arraydata, HT_F1_SD        = round(apply(arraydata[,HT_F1],1,sd), 2))
arraydata <- cbind(arraydata, GF_F1_Mean      = round(apply(arraydata[,GF_F1],1,mean), 2))
arraydata <- cbind(arraydata, GF_F1_SD        = round(apply(arraydata[,GF_F1],1,sd), 2))
arraydata <- cbind(arraydata, ratio_HT_F1_vs_GF_F1  = round(apply(arraydata[,HT_F1],1,mean) / apply(arraydata[,GF_F1],1,mean), 2))

arraydata <- cbind(arraydata, ratio_HT_BFMI_vs_HT_B6N   = round(apply(arraydata[,HT_BFMI],1,mean) / apply(arraydata[,HT_B6N],1,mean), 2))
arraydata <- cbind(arraydata, ratio_HT_BFMI_vs_HT_F1    = round(apply(arraydata[,HT_BFMI],1,mean) / apply(arraydata[,HT_F1],1,mean), 2))
arraydata <- cbind(arraydata, ratio_HT_B6N_vs_HT_F1     = round(apply(arraydata[,HT_B6N],1,mean) / apply(arraydata[,HT_F1],1,mean), 2))

arraydata <- cbind(arraydata, ratio_GF_BFMI_vs_GF_B6N   = round(apply(arraydata[,GF_BFMI],1,mean) / apply(arraydata[,GF_B6N],1,mean), 2))
arraydata <- cbind(arraydata, ratio_GF_BFMI_vs_GF_F1    = round(apply(arraydata[,GF_BFMI],1,mean) / apply(arraydata[,GF_F1],1,mean), 2))
arraydata <- cbind(arraydata, ratio_GF_B6N_vs_GF_F1     = round(apply(arraydata[,GF_B6N],1,mean) / apply(arraydata[,GF_F1],1,mean), 2))

arraydata <- cbind(arraydata, Ttest_HT_GF = apply(arraydata[,c(HT, GF)],1,function(x){
  res <- NA
  tryCatch(res <- t.test(x[1:length(HT)], x[(length(HT)+1):(length(HT)+length(GF))])$p.value,  error = function(e){ res <<- NA; })
  res
}))

arraydata <- cbind(arraydata, Ttest_HT_BFMI_GF_BFMI = apply(arraydata[,c(HT_BFMI, GF_BFMI)],1,function(x){
    res <- NA
    tryCatch(res <- t.test(x[1:length(HT_BFMI)], x[(length(HT_BFMI)+1):(length(HT_BFMI)+length(GF_BFMI))])$p.value,  error = function(e){ res <- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_HT_B6N_GF_B6N = apply(arraydata[,c(HT_B6N, GF_B6N)],1,function(x){
      res <- NA
    tryCatch(res <- t.test(x[1:length(HT_B6N)], x[(length(HT_B6N)+1):(length(HT_B6N)+length(GF_B6N))])$p.value,  error = function(e){ res <<- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_HT_B6N_HT_BFMI = apply(arraydata[,c(HT_B6N, HT_BFMI)],1,function(x){
      res <- NA
    tryCatch(res <- t.test(x[1:length(HT_B6N)], x[(length(HT_B6N)+1):(length(HT_B6N)+length(HT_BFMI))])$p.value,  error = function(e){ res <<- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_GF_B6N_GF_BFMI = apply(arraydata[,c(GF_B6N, GF_BFMI)],1,function(x){
     res <- NA
   tryCatch(res <-  t.test(x[1:length(GF_B6N)], x[(length(GF_B6N)+1):(length(GF_B6N)+length(GF_BFMI))])$p.value,  error = function(e){ res <<- NA; })
   res
}))

arraydata <- cbind(arraydata, Ttest_HT_F1_GF_F1 = apply(arraydata[,c(HT_F1, GF_F1)],1,function(x){
      res <- NA
    tryCatch(res <- t.test(x[1:length(HT_F1)], x[(length(HT_F1)+1):(length(HT_F1)+length(GF_F1))])$p.value,  error = function(e){ res <<- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_HT_F1_HT_BFMI = apply(arraydata[,c(HT_F1, HT_BFMI)],1,function(x){
     res <- NA
   tryCatch(res <-  t.test(x[1:length(HT_F1)], x[(length(HT_F1)+1):(length(HT_F1)+length(HT_BFMI))])$p.value,  error = function(e){ res <<- NA; })
   res
}))

arraydata <- cbind(arraydata, Ttest_GF_F1_GF_BFMI = apply(arraydata[,c(GF_F1, GF_BFMI)],1,function(x){
      res <- NA
    tryCatch(res <- t.test(x[1:length(GF_F1)], x[(length(GF_F1)+1):(length(GF_F1)+length(GF_BFMI))])$p.value,  error = function(e){ res <<- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_HT_F1_HT_B6N = apply(arraydata[,c(HT_F1, HT_B6N)],1,function(x){
  res <- NA
    tryCatch(res <- t.test(x[1:length(HT_F1)], x[(length(HT_F1)+1):(length(HT_F1)+length(HT_B6N))])$p.value,  error = function(e){ res <<- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_GF_F1_GF_B6N = apply(arraydata[,c(GF_F1, GF_B6N)],1,function(x){
  res <- NA
    tryCatch(res <- t.test(x[1:length(GF_F1)], x[(length(GF_F1)+1):(length(GF_F1)+length(GF_B6N))])$p.value,  error = function(e){ res <<- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_GF_F1nB6N_GF_BFMI = apply(arraydata[,c(GF_F1, GF_B6N, GF_BFMI)],1,function(x){
  res <- NA
    tryCatch(res <- t.test(x[1:(length(GF_F1)+length(GF_B6N))], x[ (length(GF_F1)+length(GF_B6N)+1):length(x) ])$p.value,  error = function(e){ res <<- NA; })
    res
}))

arraydata <- cbind(arraydata, Ttest_HT_F1nB6N_HT_BFMI = apply(arraydata[,c(HT_F1, HT_B6N, HT_BFMI)],1,function(x){
  res <- NA
    tryCatch(res <- t.test(x[1:(length(HT_F1)+length(HT_B6N))], x[ (length(HT_F1)+length(HT_B6N)+1):length(x) ])$p.value,  error = function(e){ res <<- NA; })
    res
}))

# Create some plots
#heatmap(cor(arraydata[DEtissue, arrays[,"AtlasID"]], method="spearman"))             # Correlation of samples using the probes DE between tissues
#heatmap(cor(arraydata[DEstrain, arrays[,"AtlasID"]], method="spearman"))             # Correlation of samples using the probes DE between strains

# Add the information about which probe are MultiMapping
annotationmatrix <- cbind(annotationmatrix, MultiMap = annotationmatrix[,"ProbeName"] %in% dupprobes)

if(!file.exists("Analysis/geneexpression_TestsNewNew.txt")){
  alldata <- NULL
  cnt <- 1
  ensgenes <- unique(annotationmatrix[,"ensembl_gene_id"])
  for(x in ensgenes){
    annotsubset <- annotationmatrix[which(annotationmatrix[,"ensembl_gene_id"] == x),]
    cat(paste0(cnt, "/", length(ensgenes), ", Gene:"), x, ", Probes:", dim(annotsubset)[1], "\n")
    probeinformation <- NULL
    for(x in as.character(annotsubset[,"ProbeName"])){
      probeinformation <- rbind(probeinformation, arraydata[which(rownames(arraydata) == x),])
    }
    alldata <- rbind(alldata, cbind(annotsubset, probeinformation))
    cnt <- cnt + 1
  }
  write.table(alldata, file="Analysis/geneexpression_TestsNewNew.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading gene expression data from disk\n")
  alldata <- read.table("Analysis/geneexpression_TestsNewNew.txt", sep="\t", header=TRUE)
}
