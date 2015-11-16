# Pipeline for paired end RNA-Seq analysis
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Okt, 2015
# first written Okt, 2015

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

if(!file.exists(short.gtf)){
  gtf <- read.table(reference.gtf,sep="\t")
  cat(readLines(reference.gtf,n=5),sep="\n", file=short.gtf)
  write.table(gtf[gtf[,1] %in% chrs,], file=short.gtf,sep="\t", row.names = FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
}

# Generate the command to detect SNPs in the whole population, write out as population.vcf
cat(paste0("nohup ", mpileup, " -g -f ", reference.fa," ", paste(bamfiles, collapse=" "), " -o population.bcf &"), "\n")
cat(paste0("nohup bcftools call -c population.bcf | ~/Github/bcftools/vcfutils.pl varFilter - > population.vcf &"), "\n")

# Load in the VCF data created by BCFtools
vcfdata         <- read.table("population.vcf", header = FALSE, colClasses="character")
write.table(cbind(vcfdata[, 1],vcfdata[, 2]),"SNPlocations.txt",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

vcffiles <- NULL
for(x in bamfiles){
  cat("Mpileup of file:", x, "\n")
  execute(paste0("/home/neubert/Keane/samtools-1.2/samtools mpileup -uv -t DV -t DP -l SNPlocations.txt -f ", reference.fa, " ", x, " | bcftools call -c - > ", paste0(x,".vcf")), paste0(x,".vcf"))
  vcffiles <- c(vcffiles, paste0(x,".vcf"))
}

getDP4 <- function(x){ return(sub(".*?DP4=(.*?);.*", "\\1", x)) }
if(!file.exists(paste0(readsoutput, "allsamples.dp4"))){
  dp4datafull <- NA
  for(x in vcffiles){
    cat(x, "- ")
    vcfdata <- read.table(x, sep="\t", colClasses="character")       # Read the VCF data
    vcfdata <- vcfdata[-which(grepl("INDEL", vcfdata[,8])),]         # Remove the indels
    cat("rows:", nrow(vcfdata))
    
    dp4values <- lapply(unlist(lapply(vcfdata[, 8], getDP4)), function(x){return(strsplit(x,",")[[1]])})
    ref <- unlist(lapply(dp4values,function(x){ return(as.numeric(x[1]) + as.numeric(x[2])) }))
    alt <- unlist(lapply(dp4values,function(x){ return(as.numeric(x[3]) + as.numeric(x[4])) }))
    name <- substr(strsplit(x,"/")[[1]][11], 1, 4)
    dp4 <- cbind(vcfdata[, 1], vcfdata[, 2],  vcfdata[, 4], vcfdata[, 5], ref + alt, ref, alt)
    colnames(dp4) <- c("CHROM", "POS", "REF", "ALT", paste0(name,"_DP"), paste0(name,"_REF"), paste0(name,"_ALT"))
    rownames(dp4) <- apply(dp4[,1:3],1,function(x){paste0(x,collapse="_") })
    if(is.na(dp4datafull)){                                       # First file, use all of the SNPs
      dp4datafull <- dp4
    }else{                                                        # Otherwise we take the overlap between the existing SNPs and the new SNPs
      dp4datafull <- dp4datafull[which(rownames(dp4datafull) %in% rownames(dp4)),]
      dp4 <- dp4[which(rownames(dp4) %in% rownames(dp4datafull)),]
      dp4datafull <- cbind(dp4datafull, dp4[,c(paste0(name,"_DP"), paste0(name,"_REF"), paste0(name,"_ALT"))])
      no_alt_allele <- names(which(dp4datafull[,"ALT"] == "."))
      dp4datafull[no_alt_allele,"ALT"] <- dp4[no_alt_allele,"ALT"]
    }
    cat("\nDP4:", nrow(dp4datafull),"\n")
  }
  write.table(dp4datafull, paste0(readsoutput, "allsamples.dp4"), sep="\t", quote = FALSE)
}else{
  dp4datafull <- read.table(paste0(readsoutput, "allsamples.dp4"), sep="\t", check.names=FALSE)
}

if(!file.exists(paste0(readsoutput, "allsamples_annotated.dp4"))){
  # Annotate the SNPs using the GTF data
  GTF <- read.table(short.gtf, sep="\t")                                                       # Gene models
  EXONS <- GTF[which(GTF[,3]=="exon"),]

  Gene_Exon <- t(apply(EXONS, 1, function(x){
    c(sub(".*?gene_id (.*?);.*", "\\1", x[9]), sub(".*?exon_id (.*?);.*", "\\1", x[9]), sub(".*?transcript_id (.*?);.*", "\\1", x[9]), sub(".*?exon_number (.*?);.*", "\\1", x[9]), sub(".*?gene_name (.*?);.*", "\\1", x[9]))
  }))
  Gene_Exon <- cbind(Gene_Exon, EXONS[, c(1, 4, 5)])
  colnames(Gene_Exon) <- c("gene_id", "exon_id", "transcript_id", "exon_number", "gene_name", "chr", "start", "end")

  # Annotate the SNPs using GTF information
  results <- NULL
  for(x in 1:nrow(dp4datafull)){
    snp.chr  <-  dp4datafull[x, 1]
    snp.loc  <-  as.numeric(dp4datafull[x, 2])
    cat("Going to look for a gene at:", snp.chr, snp.loc)
    idx <- which(Gene_Exon[,"chr"] == snp.chr & Gene_Exon[,"start"] <= snp.loc & Gene_Exon[,"end"] >= snp.loc)
    cat("", length(idx),"")
    if(length(idx) == 0){    # SNP is not in a gene
      results <- rbind(results, c("", "", "", "", "", "", "","", unlist(lapply(dp4datafull[x,],as.character))))
    }else if(length(idx) == 1){    # SNP is in 1 gene
      results <- rbind(results, c(unlist(lapply(Gene_Exon[idx,],as.character)), unlist(lapply(dp4datafull[x,],as.character))))
    }else{    # SNP is in multiple gene
      for(i in idx){
        results <- rbind(results, c(unlist(lapply(Gene_Exon[i,],as.character)), unlist(lapply(dp4datafull[x,],as.character))))
      }
    }
    cat("\n")
  }
  write.table(results, paste0(readsoutput, "allsamples_annotated.dp4"), sep="\t", quote = FALSE)
}else{
  results <- read.table(paste0(readsoutput, "allsamples_annotated.dp4"), sep="\t", check.names=FALSE)
}

#TODO: Fix the names of the SNPs + Add MGI ID's sice we now loose them during this step
# BIOMART annotation of the MGI_Description
if(!file.exists("/home/share/genomes/mm10/biomart/BiomartAnnotation.txt")){
  library(biomaRt)
  ensemblIDs <- as.character(dp4datafull[,"gene_id"])
  bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")                                                # Biomart for mouse genes
  biomartResults <- NULL
  for(x in seq(0, length(ensemblIDs), 1000)){                                                                        # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(ensemblIDs))                                                                       # Don't walk passed the end of the array
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

dp4datafull <- results
dp4datafull <- cbind(dp4datafull, mgi_description = NA)
for(x in 1:nrow(dp4datafull)){
  idx <- which(biomartResults[,"ensembl_gene_id"] == as.character(dp4datafull[x,"gene_id"]))
  if(length(idx) == 1) dp4datafull[x,"mgi_description"] <- as.character(biomartResults[idx, "mgi_description"])
  if(length(idx) > 1) dp4datafull[x,"mgi_description"] <- as.character(biomartResults[idx[1], "mgi_description"])
}

# Imbalance tests for the SNPs per tissue
for(x in unique(sampledescr[,"tissue"])){
  samples <- sampledescr[which(sampledescr[,"tissue"] == x), "direction"]
  names(samples) <- sampledescr[which(sampledescr[,"tissue"] == x), "Lib_id"]

  goodSNPs <- which(apply(dp4datafull[, paste0(names(samples),"_DP")],1,function(x){all(as.numeric(x) > 5)}))
  dp4dataT <- dp4datafull[goodSNPs, ]

  cat("In", x, "found", nrow(dp4dataT), "/", nrow(dp4datafull),"\n")

  B6N     <- names(samples[samples=="B6N"])
  BFMI    <- names(samples[samples=="BFMI860-12"])

  pAllele <- matrix(NA, nrow(dp4dataT), length(c(B6N,BFMI)), dimnames=list(rownames(dp4dataT), c(B6N,BFMI)))
  for(sn in c(B6N,BFMI)){
    for(snp in rownames(dp4dataT)){
      observed <- c(as.numeric(dp4dataT[snp, paste0(sn,"_REF")]), as.numeric(dp4dataT[snp, paste0(sn,"_ALT")]))
      i <- which.max(c(chisq.test(rbind(observed, c(sum(observed), 1)))$p.value, chisq.test(observed)$p.value, chisq.test(rbind(observed, c(1, sum(observed))))$p.value))
      pAllele[snp, sn] <- c("Ref", "Hetro", "Alt")[i]
    }
  }

  cat("Parental alleles determined\n")

  ratios  <- matrix(NA, nrow(dp4dataT), length(samples), dimnames=list(rownames(dp4dataT),names(samples)))
  probs   <- matrix(NA, nrow(dp4dataT), length(samples), dimnames=list(rownames(dp4dataT),names(samples)))
  chisqs  <- matrix(NA, nrow(dp4dataT), length(samples), dimnames=list(rownames(dp4dataT),names(samples)))
  chisqsS <- matrix(NA, nrow(dp4dataT), length(samples), dimnames=list(rownames(dp4dataT),names(samples)))
  for(sn in names(samples)){
    ratios[,sn]  <- round((as.numeric(dp4dataT[, paste0(sn,"_ALT")]) / as.numeric(dp4dataT[, paste0(sn,"_DP")])) * 100,d=1)
    chiTests     <- apply(dp4dataT[, c(paste0(sn,"_REF"), paste0(sn,"_ALT"))], 1, function(observed){chisq.test(as.numeric(observed))})
    probs[,sn]   <- unlist(lapply(chiTests, function(x){x$p.value}))
    chisqs[,sn]  <- round(sign(ratios[,sn] - 50)  * as.numeric(unlist(lapply(chiTests, function(x){unlist(x)["statistic.X-squared"]}))), d = 2)
    chisqsS[,sn] <- round(sign(ratios[,sn] - 50)  * sqrt(as.numeric(unlist(lapply(chiTests, function(x){unlist(x)["statistic.X-squared"]})))), d = 2)
  }
  colnames(pAllele) <- paste0(samples[colnames(pAllele)], "_", colnames(pAllele),"_call")
  colnames(ratios)  <- paste0(samples[colnames(ratios)], "_", colnames(ratios))
  colnames(probs)   <- paste0(samples[colnames(probs)], "_", colnames(probs), "_p")
  colnames(chisqs)  <- paste0(samples[colnames(chisqs)], "_", colnames(chisqs), "_chisq")
  colnames(chisqsS)  <- paste0(samples[colnames(chisqsS)], "_", colnames(chisqsS), "_chisqS")

  cat("Allelic imbalanced detected\n")
  output <- cbind(dp4dataT[,9:12], dp4dataT[,1:5], dp4dataT[,"mgi_description"], dp4dataT[,6:8], pAllele, ratios, probs, chisqs, chisqsS)
  write.table(output, file=paste0(readsoutput, "snp_stats_",x,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}

for(x in unique(sampledescr[,"tissue"])){
  samples <- sampledescr[which(sampledescr[,"tissue"] == x), "direction"]
  names(samples) <- sampledescr[which(sampledescr[,"tissue"] == x), "Lib_id"]
  
  B6N     <- names(samples[samples=="B6N"])
  BFMI    <- names(samples[samples=="BFMI860-12"])
  matB6N  <- names(samples[samples=="B6NxBFMI860-12"])
  matBFMI <- names(samples[samples=="BFMI860-12xB6N"])
  
  snpstats <- read.csv(paste0(readsoutput, "snp_stats_",x,".txt"), sep="\t", check.names=FALSE, header=TRUE)
  calldata <- snpstats[,paste0(samples[c(B6N,BFMI)], "_", c(B6N,BFMI),"_call")]
  pDiff <- which(apply(calldata,1,function(x){                                        # Parents need to be homozygous and different
    if(any(x == "Hetro")) return(FALSE)
    if(x[1] == x[2] && x[3] == x[4] && x[1] != x[3]) return(TRUE)
    return(FALSE)
  }))
  cat("Different between parentals", nrow(snpstats[pDiff,]), "/", nrow(snpstats), "\n")

  snpstats <- snpstats[pDiff,]

  matB6N.mean  <- apply(snpstats[, paste0(samples[matB6N], "_", matB6N)], 1, function(x){return(mean(x))})
  matB6N.sd    <- apply(snpstats[, paste0(samples[matB6N], "_", matB6N)], 1, function(x){return(sd(x))})
  matBFMI.mean <- apply(snpstats[, paste0(samples[matBFMI], "_", matBFMI)], 1, function(x){return(mean(x))})
  matBFMI.sd   <- apply(snpstats[, paste0(samples[matBFMI], "_", matBFMI)], 1, function(x){return(sd(x))})

  
  matB6N.sum  <- apply(snpstats[, paste0(samples[matB6N], "_", matB6N, "_chisqS")], 1, function(x){return(sum(x))})
  matBFMI.sum <- apply(snpstats[, paste0(samples[matBFMI], "_", matBFMI, "_chisqS")], 1, function(x){return(sum(x))})
  chisq.diff  <- apply(cbind(matBFMI.sum,matB6N.sum),1,diff)

  snpstats <- cbind(snpstats, matBFMI.mean, matBFMI.sd, matB6N.mean, matB6N.sd, matBFMI.sum, matB6N.sum, chisq.diff)
  
  idx <- which(snpstats[,"matB6N.sd"] > 20 | snpstats[,"matBFMI.sd"] > 20)
  cat("Removing", length(idx), "SNPs which vary too much ( > 20 %)\n")
  snpstats <- snpstats[-idx,]

  ## Permutation to find which delta Chi2 scores are significantly different
  randomScoresSum <- NULL ; randomScoresDiff <- NULL
  for(l in 1:10000){
    i <- sample(nrow(snpstats), 1)
    randomScoresSum <- c(randomScoresSum, max(abs(c(snpstats[i,"matB6N.sum"], snpstats[i,"matBFMI.sum"])),na.rm=TRUE))
    randomScoresDiff <- c(randomScoresDiff, max(abs(snpstats[i,"chisq.diff"]),na.rm=TRUE))
  }
  threshold.a5s <- sort(randomScoresSum)[length(randomScoresSum) * .95] ; threshold.a5d <- sort(randomScoresDiff)[length(randomScoresDiff) * .95]

  cat(x, " - Theshold single",  threshold.a5s, ", Threshold Difference", threshold.a5d , "\n")
  
  write.table(snpstats[which(snpstats[,"matB6N.sum"] > threshold.a5s | snpstats[,"matB6N.sum"] < -threshold.a5s), ], file=paste0(readsoutput, "snp_stats_",x,"_ASE_mB6N.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snpstats[which(snpstats[,"matBFMI.sum"] > threshold.a5s | snpstats[,"matBFMI.sum"] < -threshold.a5s), ], file=paste0(readsoutput, "snp_stats_",x,"_ASE_mBFMI.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(snpstats[which(snpstats[,"chisq.diff"] > threshold.a5d | snpstats[,"chisq.diff"] < -threshold.a5d),], file=paste0(readsoutput, "snp_stats_",x,"_ASE_diff.txt"), sep="\t", quote=FALSE, row.names=FALSE)
}

