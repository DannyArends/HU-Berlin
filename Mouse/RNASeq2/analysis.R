# Pipeline for paired end RNA-Seq analysis
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Okt, 2015
# first written Okt, 2015

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

baminputdir     <- as.character("/home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/FASTQ/")
readsoutput     <- as.character("/home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/Analysis/")
reference.gtf   <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.81.gtf"                      # wget ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz
short.gtf       <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.81.short.gtf"
genesonly.gtf   <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.81.genesonly.gtf"
sampledescrpath <- "/home/arends/RNASeq/Experiment/SampleDescription.txt"
chrdescrpath    <- "/home/arends/Mouse/GTF/MouseChrInfo.txt"
reference       <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.dna"
reference.fa    <- paste0(reference, ".fa")        #!
mpileup         <- "~/Github/samtools/samtools mpileup"

chrs <- c(1:19,"X","Y","MT")
sampledescr  <- read.table(sampledescrpath, sep="\t", header=TRUE, colClasses="character")
chrominfo    <- read.table(chrdescrpath, sep="\t", header=TRUE, colClasses=c("character","integer","logical"))
chrominfo    <- chrominfo[chrominfo[,1] %in% chrs,]

baminputfolder <- dir(baminputdir)
inputbfiles <- baminputfolder[which(grepl("output", baminputfolder))]
sampleNames <- substr(inputbfiles, 1, 4)
cat("Found ", length(sampleNames), " samples in folder: ", baminputdir, "\n")

bamfiles <- NULL
for(msample in inputbfiles){
  files <- dir(paste0(baminputdir, msample))
  fname <- files[which(grepl(".aln.sort.rmdup.rg.realigned.recal.bam$", files))]
  if(!((length(fname) == 0) && (typeof(fname) == "character"))) bamfiles <- c(bamfiles, paste0(baminputdir, msample, "/", fname))
}                                                                                                             #done <- substr(bamfiles, 1, 4)
cat("Found ", length(bamfiles), " aligned bam files in folder: ", baminputdir, "\n")

completed <- unlist(lapply(strsplit(unlist(lapply(strsplit(bamfiles,"/"),"[",10)),"_"),"[",1))


if(!file.exists(short.gtf) || if(!file.exists(genesonly.gtf))) gtf <- read.table(reference.gtf,sep="\t")
if(!file.exists(short.gtf)){
  cat(readLines(reference.gtf,n=5),sep="\n", file=short.gtf)
  write.table(gtf[gtf[,1] %in% chrs,], file=short.gtf,sep="\t", row.names = FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
}
if(!file.exists(genesonly.gtf)){  
  cat(readLines(reference.gtf,n=5),sep="\n", file=genesonly.gtf)
  okG <- grepl("ensembl", gtf[,2]) | grepl("insdc", gtf[,2]) | grepl("mirbase", gtf[,2])
  write.table(gtf[gtf[,1] %in% chrs & gtf[,3] == "gene" & okG,], file=genesonly.gtf, sep="\t", row.names = FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
}

# # Use gtf2bed to create the reference transcriptome in BED format
reference.bed <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.81.short.bed"
execute(paste0("PATH=/home/arends/HLRN/bin:$PATH && gtf2bed < ", short.gtf, " > ", reference.bed), reference.bed)

# # Use gtf2bed to create the reference transcriptome in BED format
genesonly.bed <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.81.genesonly.bed"
execute(paste0("PATH=/home/arends/HLRN/bin:$PATH && gtf2bed < ", genesonly.gtf, " > ", genesonly.bed), genesonly.bed)


# Create index file symbolic links, to adhere to the bedtools expectation
#for(x in bamfiles){
#  old <- gsub("bam", "bai", x); new <- gsub("bam", "bam.bai", x);  out <- gsub("bam", "counts", x)
#  execute(paste0("ln -s ",old," ", new, "\n"), new)
#  execute(paste0("~/HLRN/bedtools2/bin/bedtools multicov -bams ", x, " -bed /home/share/genomes/mm10/Mus_musculus.GRCm38.81.short.bed > ", out), out)
#}

# Extract the counts from the individual bam files
out <- paste0(readsoutput, "raw_counts_genesonly.txt")
execute(paste0("nohup ~/HLRN/bedtools2/bin/bedtools multicov -split -q 10 -bams ", paste0(bamfiles, collapse=" "), " -bed /home/share/genomes/mm10/Mus_musculus.GRCm38.81.genesonly.bed > ", out, " &"), out)

if(!file.exists(paste0(readsoutput, "RPKM_norm_log.txt"))){
  expressiondata <- read.table(out, sep="\t")
  expressiondata <- expressiondata[,-c(4,5,7,8,9)]
  geneinfo <- unlist(lapply(expressiondata[,5],as.character))
  ensmusg <- gsub("gene_id ","",unlist(lapply(strsplit(geneinfo,"; "),"[",1)))
  gname <- gsub("gene_name ","",unlist(lapply(strsplit(geneinfo,"; "),"[",3)))

  samplenames <- unlist(lapply(strsplit(unlist(lapply(strsplit(bamfiles,"/"),"[",10)),"_"),"[",1))
  expressiondata <- cbind(expressiondata[,c(1,2,3,4)], ensmusg, gname, expressiondata[,c(6:ncol(expressiondata))])
  colnames(expressiondata) <- c("Chr", "Start", "End", "Strand", "ensembl_gene_id", "gene_name", samplenames)

  rawreads <- expressiondata[,samplenames]

  rawreads[rawreads == 0] <- NA                                                               # Change 0 reads to NA, so we can do quantile normalisation
  rawreadsQnorm <- normalize.quantiles(as.matrix(rawreads))
  rawreadsQnorm[is.na(rawreadsQnorm)] <- 0                                                    # Change back the NAs to 0 reads

  # RPKM = (10^9 * C)/(N * L)
  # C = Number of reads mapped to a gene
  # N = Total mapped reads in the sample
  # L = gene length in base-pairs for a gene

  chrominfo     <- read.table("/home/share/genomes/mm10/MouseChrInfo.txt", sep="\t", header=TRUE, colClasses=c("character","integer","logical"))
  mouse         <- makeTranscriptDbFromGFF(short.gtf, format = "gtf", exonRankAttributeName="exon_number", 
                                           species="Mus musculus", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/")

                                           exonsByGene   <- exonsBy(mouse, by = "gene")
  exonicGeneSizes <- lapply(exonsByGene, function(x){ sum(width(reduce(x))) })                # Get the length of each gene using only the exons
  N <- apply(rawreadsQnorm, 2, sum)                                                           # Get the number of reads in all samples

  orderedSizes <- exonicGeneSizes[as.character(expressiondata[,"ensembl_gene_id"])]

  n <- 1
  RPKM <- t(apply(rawreadsQnorm, 1, function(C){                                              # Calculate the RPKM values per gene
    L     <- as.numeric(orderedSizes[n])
    RPKM  <- (10^9 * C) / (N * L)
    n    <<- n + 1
    return(round(RPKM, d = 3))
  }))
  RPKM <- apply(RPKM, 2, function(x){log(x+1)})
  RPKM <- apply(RPKM, 2, function(x){round(x,3)})
  expressiondata[,samplenames] <- RPKM

  cat("We called expressions for", nrow(expressiondata), "genes\n")
  write.table(expressiondata, file=paste0(readsoutput, "RPKM_norm_log.txt"), sep="\t",row.names=FALSE)                                  # Write the normalized RPKM values to file
}else{
  expressiondata <- read.table(paste0(readsoutput, "RPKM_norm_log.txt"), sep="\t",header=TRUE,check.names=FALSE)
}

# Calculate statistics
statistics <- NULL
for(x in unique(sampledescr[,"tissue"])){
  samples <- sampledescr[which(sampledescr[,"tissue"] == x), "direction"]
  names(samples) <- sampledescr[which(sampledescr[,"tissue"] == x), "Lib_id"]
  #expressed <- which(!apply(expressiondata[, names(samples)], 1, function(x){all(as.numeric(x) < 0.5)}))
  #cat("Tissue:",x, "has", length(expressed),"expressed\n")
  #expSubset <- expressiondata[expressed,]
  
  B6N     <- names(samples[samples=="B6N"])
  BFMI    <- names(samples[samples=="BFMI860-12"])
  matB6N  <- names(samples[samples=="B6NxBFMI860-12"])
  matBFMI <- names(samples[samples=="BFMI860-12xB6N"])
  
  meanB6N      <- round(apply(expressiondata[, B6N], 1, function(x){mean(x)}), d = 2)
  meanBFMI     <- round(apply(expressiondata[, BFMI], 1, function(x){mean(x)}), d = 2)
  meanMatB6N   <- round(apply(expressiondata[, matB6N], 1, function(x){mean(x)}), d = 2)
  meanMatBFMI  <- round(apply(expressiondata[, matBFMI], 1, function(x){mean(x)}), d = 2)
  pvalues      <- apply(expressiondata[, c(matB6N,matBFMI)], 1, function(x){
    tryCatch(return(t.test(x[1:3],x[4:6])$p.value),  error = function(e){ return(NA) } )
  })
  ratios       <- round(meanMatBFMI / meanMatB6N, d=2)
  
  expdata <- cbind(meanB6N, meanBFMI, meanMatB6N, meanMatBFMI, pvalues, ratios)
  colnames(expdata) <- paste0(c("mean_B6N", "mean_BFMI", "mean_B6NxBFMI860-12", "mean_BFMI860-12xB6N", "pValues", "ratios"), "_", substr(x,1,1))
  statistics <- cbind(statistics, expdata)
}

allsamples <- paste0(sampledescr[,"direction"], "_", substr(sampledescr[,"tissue"],1,1))
names(allsamples) <- sampledescr[, "Lib_id"]

idx    <- which(sampledescr[,"direction"] %in% c("BFMI860-12xB6N","B6NxBFMI860-12"))
snames <- sampledescr[idx, "Lib_id"]
anovas <- t(apply(expressiondata[,snames],1,function(ge){
  return(anova(lm(ge ~ sampledescr[idx,"tissue"] + sampledescr[idx,"direction"]))[[5]][1:2])
}))
colnames(anovas) <- c("pValues_Tissue", "pValues_DOC")

colnames(expressiondata)[7:ncol(expressiondata)] <- allsamples[colnames(expressiondata)[7:ncol(expressiondata)]]

# BIOMART annotation of the MGI_Description
if(!file.exists("/home/share/genomes/mm10/biomart/BiomartAnnotation.txt")){
  library(biomaRt)
  ensemblIDs <- as.character(expressiondata[,"ensembl_gene_id"])
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")                                                # Biomart for mouse genes
  biomartResults <- NULL
  for(x in seq(0, length(ensemblIDs), 1000)){                                                                     # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(ensemblIDs))                                                                    # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")

    res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("ensembl_gene_id"), 
                        values = ensemblIDs[x:xend], mart = bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
    Sys.sleep(1)
  }
  write.table(biomartResults, file="/home/share/genomes/mm10/biomart/BiomartAnnotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading biomart annotation from disk\n")
  biomartResults <- read.table("/home/share/genomes/mm10/biomart/BiomartAnnotation.txt", sep="\t", header=TRUE, colClasses=c("character"))
}

# Add MGI description to the expression data
expressiondata <- cbind(expressiondata[,1:6], mgi_description = NA, expressiondata[,7:ncol(expressiondata)])
for(x in 1:nrow(expressiondata)){
  idx <- which(biomartResults[,"ensembl_gene_id"] == expressiondata[x,"ensembl_gene_id"])
  if(length(idx) == 1) expressiondata[x,"mgi_description"] <- biomartResults[idx, "mgi_description"]
  if(length(idx) > 1) expressiondata[x,"mgi_description"] <- biomartResults[idx[1], "mgi_description"]
}

# Add statistics to the expression data
expdatastats <- cbind(expressiondata, statistics, anovas)

# SNPs detected by DNA-Seq
bfmi.snps <- "/home/share/genomes/mm10/SNPs_BFMI/860v2.pes.s.rmd.rg.realigned.recal.haplotclld.Q30DP3HomFilt_NovelSNPer_detailed.out"
snps <- read.csv(bfmi.snps, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors = FALSE)   # DNA-Seq SNP data
snps <- snps[-which(snps$FuncClass == "SYNONYMOUS_CODING"), ]                                     # Throw away the SYNONYMOUS_CODING functional class

codingRegion  <- unique(snps[grepl("CODING_REGION", snps$Region),"Gene_ID"])
domainRegion  <- unique(snps[!is.na(snps$DomainRegion),"Gene_ID"])
UTR5          <- unique(snps[grepl("5PRIME_UTR", snps$Region),"Gene_ID"])
UTR3          <- unique(snps[grepl("3PRIME_UTR", snps$Region),"Gene_ID"])
SPLICE        <- unique(snps[grepl("SPLICE_SITE", snps$Region),"Gene_ID"])
STARTSTOP     <- unique(snps[grepl("START", snps$FuncClass),"Gene_ID"])
STARTSTOP     <- unique(c(STARTSTOP, unique(snps[grepl("STOP", snps$FuncClass),"Gene_ID"])))

expdatastats <- cbind(expdatastats, "CodingRegion"  = rep(0, nrow(expdatastats))); expdatastats[which(expdatastats$ensembl_gene_id %in% codingRegion), "CodingRegion"] <- 1
expdatastats <- cbind(expdatastats, "DomainRegion"  = rep(0, nrow(expdatastats))); expdatastats[which(expdatastats$ensembl_gene_id %in% domainRegion), "DomainRegion"] <- 1
expdatastats <- cbind(expdatastats, "5' UTR"        = rep(0, nrow(expdatastats))); expdatastats[which(expdatastats$ensembl_gene_id %in% UTR5), "5' UTR"] <- 1
expdatastats <- cbind(expdatastats, "3' UTR"        = rep(0, nrow(expdatastats))); expdatastats[which(expdatastats$ensembl_gene_id %in% UTR3), "3' UTR"] <- 1
expdatastats <- cbind(expdatastats, "SpliceSite"    = rep(0, nrow(expdatastats))); expdatastats[which(expdatastats$ensembl_gene_id %in% SPLICE), "SpliceSite"] <- 1
expdatastats <- cbind(expdatastats, "Stop/Start"    = rep(0, nrow(expdatastats))); expdatastats[which(expdatastats$ensembl_gene_id %in% STARTSTOP), "Stop/Start"] <- 1

write.table(expdatastats, file=paste0(readsoutput, "RPKM_norm_log_stats.txt"), sep="\t",row.names=FALSE)                                  # Write the normalized RPKM values to file


