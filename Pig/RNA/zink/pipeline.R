# Pipeline for paired end RNA-Seq analysis
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2015
# first written Aug, 2015

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping this step\n"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish\n\n"); q("no") }
}

# Input filename (Should be supplied as input to this script)          e.g.:       nohup Rscript pipeline.R /home/arends/RNASeq/FastQ/4422_GCCAAT_L001 > 4422_GCCAAT_L001.log 2>&1&
# Input filename (Should be supplied as input to this script)          e.g.:       nohup Rscript pipeline.R /home/arends/NAS/Mouse/RNA/Sequencing/ReciprocalCrossB6xBFMI/FASTQ/4422_GCCAAT_L001 > 4422_GCCAAT_L001.log 2>&1&
inputfile <- as.character(commandArgs(trailingOnly = TRUE)[1])
#inputfile <-"/home/arends/RNASeq/FastQ/4422_GCCAAT_L001"

# Create an output folder, and a store this as the output filebase
dir.create(paste0(inputfile,".output/"))
base         <- strsplit(inputfile, "/")[[1]]
fname        <- base[length(base)]
filebase     <- paste0(paste0(inputfile,".output/"), fname)
readgroupID  <- as.numeric(fname)

cat("base:", base, "\n"); cat("fname:", fname, "\n"); cat("filebase:", filebase, "\n"); cat("readgroupID:", readgroupID, "\n")

reference              <- "/home/share/genomes/Sscrofa10.2/Sus_scrofa.Sscrofa10.2.dna"
reference.fa           <- paste0(reference, ".fa")        #!
reference.fa.gz        <- paste0(reference, ".fa.gz")     #!
reference.fa.fai       <- paste0(reference, ".fa.fai")
reference.dict         <- paste0(reference, ".dict")
reference.bt2idx       <- paste0(reference, ".bt2idx")
reference.bt2idx.fa    <- paste0(reference, ".bt2idx.fa") #!

# wget ftp://ftp.ensembl.org/pub/release-84/variation/vcf/sus_scrofa/Sus_scrofa.vcf.gz
# Recode using VCF tools to remove info fields (containing white spaces)
# ./vcftools --gzvcf /home/share/genomes/Sscrofa10.2/Sus_scrofa.vcf.gz --remove-indels --recode -c | bgzip -c > /home/share/genomes/Sscrofa10.2/Sus_scrofa.recode.vcf.gz

reference.snps      <- "/home/share/genomes/Sscrofa10.2/Sus_scrofa.recode.vcf.gz"
reference.gtf       <- "/home/share/genomes/Sscrofa10.2/Sus_scrofa.Sscrofa10.2.84.gtf"                   # wget ftp://ftp.ensembl.org/pub/release-84/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.84.gtf.gz

picarddict.exec     <- "java -Xmx4g -jar /opt/picard-tools-1.99/CreateSequenceDictionary.jar"
picardrmdup.exec    <- "java -Xmx4g -jar /opt/picard-tools-1.99/MarkDuplicates.jar"
picardrg.exec       <- "java -Xmx4g -jar /opt/picard-tools-1.99/AddOrReplaceReadGroups.jar"
picardreo.exec      <- "java -Xmx4g -jar /opt/picard-tools-1.99/ReorderSam.jar"
gatk.exec           <- "java -Xmx16g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"
b2build.exec        <- "bowtie2-build"
bowtie2.exec        <- "bowtie2"
samtools.exec       <- "samtools"
tabix.exec          <- "tabix"


# Download only the chromosomes we want from the reference (if it does not exist)
if(!file.exists(reference.fa)){
  chromosomes <- c(1:18, "X", "Y", "MT")
  for(chr in chromosomes){ execute(paste0("curl ftp://ftp.ensembl.org/pub/release-84/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa10.2.dna.chromosome.", chr, ".fa.gz | gunzip -c >> ", reference.fa)) }
}else{ cat("Reference was found at", reference.fa,"\n") }

if(!file.exists(reference.fa.gz)) execute(paste0("gzip < ", reference.fa, " > ", reference.fa.gz), reference.fa.gz)

# Index the reference genome
execute(paste0(picarddict.exec, " R=", reference.fa, " O=", reference.dict), reference.dict)
execute(paste0(b2build.exec, " -f ", reference.fa, " ", reference.bt2idx), paste0(reference.bt2idx,".1.bt2"))
execute(paste0(samtools.exec, " faidx ", reference.fa), reference.fa.fai)

#execute(paste0(tabix.exec, " ", reference.indels), paste0(reference.indels,".tbi"))
execute(paste0(tabix.exec," ", reference.snps), paste0(reference.snps,".tbi"))

# For tophat2 create a symbolic link to the reference fasta file
execute(paste0("ln -s ",  reference.fa, " ", reference.bt2idx.fa), reference.bt2idx.fa)

# Create tophat2 index file for gtf
transcriptome.index  <- "/home/share/genomes/Sscrofa10.2/tophat2/Sscrofa10.2.genes"
execute(paste0("tophat2 -G ", reference.gtf, " --transcriptome-index=", transcriptome.index, " ", reference.bt2idx), transcriptome.index)

# Trim reads by Trimmomatic
trimmomatic.files  <- c(paste0(inputfile,"_R1.fastq.gz"), paste0(inputfile,"_R2.fastq.gz"), paste0(filebase,"_R1.P.fastq.gz"), paste0(filebase,"_R1.U.fastq.gz"), paste0(filebase,"_R2.P.fastq.gz"), paste0(filebase,"_R2.U.fastq.gz"))
trimmomatic.exec   <- "java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar"
trimmomatic.opts   <- "ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

execute(paste0(trimmomatic.exec, " PE ", paste0(trimmomatic.files, collapse=" "), " ", trimmomatic.opts), trimmomatic.files[3])

# Alignment using bowtie2
#alignment.file <- paste0(filebase, ".aln.bam")
#execute(paste0(bowtie2.exec, " -x ", reference.bt2idx, " -1 ", trimmomatic.files[3]," -2 ", trimmomatic.files[5]," -U ", trimmomatic.files[4], ",", trimmomatic.files[6]," -X 2000 -I 50 | ",samtools.exec," view -bS - > ",alignment.file), alignment.file)

# Aligment using tophat2
alignment.file    <- paste0(filebase, ".tophat2.aln.bam")
tophat.options  <- paste0("--prefilter-multihits --rg-id ", readgroupID, " --rg-sample ", fname, " --rg-library RNA-seq --rg-platform Illumina --b2-sensitive --read-gap-length 12 --read-mismatches 2 --read-edit-dist 13 --max-insertion-length 20 --max-deletion-length 30 --num-threads 2",collapse="")
tophat.log      <- paste0(filebase, ".tophat2.log")
execute(paste0("tophat2"," -o ", alignment.file, " --transcriptome-index ", transcriptome.index, " ", tophat.options, " ", reference.bt2idx, " ", trimmomatic.files[3], " ", trimmomatic.files[5], ",", trimmomatic.files[4], ",", trimmomatic.files[6], " > ", tophat.log, collapse=""), alignment.file)
alignment.file    <- paste0(filebase, ".tophat2.aln.bam/accepted_hits.bam")

# Sort and remove duplicates
samtools.file    <- paste0(filebase, ".aln.sort")
samtools.bamfile <- paste0(samtools.file, ".bam")
picard.file <- paste0(filebase, ".aln.sort.rmdup.bam")
metrics.file <- paste0(filebase, ".metrics.txt")

execute(paste0(samtools.exec, " sort -@ 4 -m 2G ",alignment.file, " ", samtools.file), samtools.bamfile)
execute(paste0(picardrmdup.exec, " REMOVE_DUPLICATES=true INPUT=", samtools.bamfile, " OUTPUT=", picard.file, " METRICS_FILE=", metrics.file), picard.file)
#execute(paste0(samtools.exec, " flagstat ", picard.file))

# Add read groups
picard.rgfile             <-  paste0(filebase, ".aln.sort.rmdup.rg.bam")
picard.rgfile.reordered   <-  paste0(filebase, ".aln.sort.rmdup.rg.o.bam")
execute(paste0(picardrg.exec, " INPUT=", picard.file, " OUTPUT=", picard.rgfile, " CREATE_INDEX=false RGID=",readgroupID," RGLB=LIB860 RGPL=Illumina RGPU=X RGSM=860"), picard.rgfile)

# index the bamfile
indexed.rgfile <-  paste0(filebase, ".aln.sort.rmdup.rg.bam.bai")
execute(paste0(samtools.exec, " index ", picard.rgfile), indexed.rgfile)

# Reorder the bam file so it matches the ordering of the fasta file (how it gets to be out of order, I have no f*cking clue)
execute(paste0(picardreo.exec, " I=", picard.rgfile," O=", picard.rgfile.reordered," REFERENCE=",reference.fa), picard.rgfile.reordered)

# index the bamfile
indexed.rgfile <-  paste0(filebase, ".aln.sort.rmdup.rg.o.bam.bai")
execute(paste0(samtools.exec, " index ", picard.rgfile.reordered), indexed.rgfile)

# Base Recalibration
gatk.cov1 <- paste0(filebase, ".1.covariates")
execute(paste0(gatk.exec, " -nct 4 -T BaseRecalibrator -R ", reference.fa," -knownSites ", reference.snps, " -I ", picard.rgfile.reordered, " -o ", gatk.cov1," -U ALLOW_N_CIGAR_READS"), gatk.cov1)
gatk.recal <- paste0(filebase, ".aln.sort.rmdup.rg.realigned.recal.bam")
execute(paste0(gatk.exec, " -nct 4 -T PrintReads -R ", reference.fa," -I ", picard.rgfile.reordered," -BQSR ", gatk.cov1," -U ALLOW_N_CIGAR_READS -o ", gatk.recal), gatk.recal)
gatk.cov2 <- paste0(filebase, ".2.covariates")
execute(paste0(gatk.exec, " -nct 4 -T BaseRecalibrator -R ", reference.fa," -knownSites ", reference.snps, "  -I ",gatk.recal," -o ", gatk.cov2," -U ALLOW_N_CIGAR_READS"), gatk.cov2)
gatk.plots <- paste0(filebase, ".plots.pdf")
execute(paste0(gatk.exec, " -T AnalyzeCovariates -R ", reference.fa," -before ", gatk.cov1, " -after ", gatk.cov2," -U ALLOW_N_CIGAR_READS -plots ", gatk.plots), gatk.plots)

q("no")
