# Pipeline for paired end RNA-Seq analysis
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2015
# first written Aug, 2015

# Execute function, does not execute when outputfile exists
execute <- function(x, outputfile = NA, intern = FALSE){
  if(!is.na(outputfile) && file.exists(outputfile)){ cat("Output for step exists, skipping"); return("") }
  cat("----", x, "\n"); res <- system(x, intern = intern); cat(">>>>", res[1], "\n")
  if(res[1] >= 1){ cat("Error external process did not finish"); q("no") }
}

# Input filename (Should be supplied as input to this script)          e.g.:       nohup Rscript pipeline.R /home/arends/RNASeq/FastQ/4422_GCCAAT_L001 > 4422_GCCAAT_L001.log 2>&1&
inputfile <- as.character(commandArgs(trailingOnly = TRUE)[1])

# Create an output folder, and a store this as the output filebase
dir.create(paste0(inputfile,".output/"))
base       <- strsplit(inputfile, "/")[[1]]
filebase   <- paste0(paste0(inputfile,".output/"), base[length(base)])

# Locations, executable names and constants
sampledescription   <- "/home/arends/RNASeq/Experiment/SampleDescription.txt"
reference           <- "/home/share/genomes/mm10/Mus_musculus.GRCm38.dna"
reference.fa        <- paste0(reference, ".fa")
reference.fa.gz     <- paste0(reference, ".fa.gz")
reference.fa.fai    <- paste0(reference, ".fa.fai")
reference.dict      <- paste0(reference, ".dict")
reference.bt2idx    <- paste0(reference, ".bt2idx")

reference.snps      <- "/home/share/genomes/mm10/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz"
reference.indels    <- "/home/share/genomes/mm10/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz"
indels.intervals    <- "/home/share/genomes/mm10/output.indels.intervals"

picarddict.exec     <- "java -Xmx4g -jar /opt/picard-tools-1.99/CreateSequenceDictionary.jar"
picardrmdup.exec    <- "java -Xmx4g -jar /opt/picard-tools-1.99/MarkDuplicates.jar"
picardrg.exec       <- "java -Xmx4g -jar /opt/picard-tools-1.99/AddOrReplaceReadGroups.jar"
gatk.exec           <- "java -Xmx4g -jar /opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"
b2build.exec        <- "bowtie2-build"
bowtie2.exec        <- "bowtie2"
samtools.exec       <- "samtools"
tabix.exec          <- "tabix"

# Download only the chromosomes we want from the reference (if it does not exist)
if(!file.exists(reference.fa)){
  chromosomes <- c(1:19, "X", "Y", "MT")
  for(chr in chromosomes){ execute(paste0("curl ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.", chr, ".fa.gz | gunzip -c >> ", reference.fa)) }
}else{ cat("Reference was found at", reference.fa,"\n") }
# Create a gzip version of the reference
execute(paste0("gzip < ", reference.fa, " > ", reference.fa.gz), reference.fa.gz)

# Index the reference genome
execute(paste0(picarddict.exec, "R=",reference.fa.gz," O=", reference.dict), reference.dict)
execute(paste0(b2build.exec, " -f ", reference.fa.gz, " ", reference.bt2idx), paste0(reference.bt2idx,".1.bt2"))
execute(paste0(samtools.exec, " faidx ", reference.fa.gz), reference.fa.fai)

execute(paste0(tabix.exec, " ", reference.indels), paste0(reference.indels,".tbi"))
execute(paste0(tabix.exec," ", reference.snps), paste0(reference.snps,".tbi"))

# Create indel realignment file
execute(paste0(gatk.exec, " -nt 4 -T RealignerTargetCreator -R ", reference.fa, " -known ", reference.indels," -o ", indels.intervals, " -U ALLOW_N_CIGAR_READS"), indels.intervals)

# Trimmomatic
trimmomatic.files  <- c(paste0(inputfile,"_R1_001.subset.fastq.gz"), paste0(inputfile,"_R2_001.subset.fastq.gz"), paste0(filebase,"_R1.P.fastq.gz"), paste0(filebase,"_R1.U.fastq.gz"), paste0(filebase,"_R2.P.fastq.gz"), paste0(filebase,"_R2.U.fastq.gz"))
trimmomatic.exec   <- "java -jar /opt/Trimmomatic-0.32/trimmomatic-0.32.jar"
trimmomatic.opts   <- "ILLUMINACLIP:/opt/Trimmomatic-0.32/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

execute(paste0(trimmomatic.exec, " PE ", paste0(trimmomatic.files, collapse=" "), " ", trimmomatic.opts), trimmomatic.files[3])

# Alignment
bowtie2.file <- paste0(filebase, ".aln.bam")
execute(paste0(bowtie2.exec, " -x ", reference.bt2idx, " -1 ", trimmomatic.files[3]," -2 ", trimmomatic.files[5]," -U ", trimmomatic.files[4], ",", trimmomatic.files[6]," -X 2000 -I 50 | ",samtools.exec," view -bS - > ",bowtie2.file), bowtie2.file)

# Sort and remove duplicates
samtools.file    <- paste0(filebase, ".aln.sort")
samtools.bamfile <- paste0(samtools.file, ".bam")
picard.file <- paste0(filebase, ".aln.sort.rmdup.bam")
metrics.file <- paste0(filebase, ".metrics.txt")

execute(paste0(samtools.exec, " sort -@ 4 -m 2G ",bowtie2.file, " ", samtools.file), samtools.bamfile)
execute(paste0(picardrmdup.exec, " REMOVE_DUPLICATES=true INPUT=", samtools.bamfile, " OUTPUT=", picard.file, " METRICS_FILE=", metrics.file), picard.file)
execute(paste0(samtools.exec, " flagstat ", picard.file))

# Add read groups
picard.rgfile <-  paste0(filebase, ".aln.sort.rmdup.rg.bam")
execute(paste0(picardrg.exec, " INPUT=", picard.file, " OUTPUT=", picard.rgfile, " CREATE_INDEX=false RGID=860 RGLB=LIB860 RGPL=Illumina RGPU=X RGSM=860"), picard.rgfile)

# index the bamfile
indexed.rgfile <-  paste0(filebase, ".aln.sort.rmdup.rg.bam.bai")
execute(paste0(samtools.exec, " index ", picard.rgfile), indexed.rgfile)

# Indel realign
gatk.realigned <- paste0(filebase, ".aln.sort.rmdup.rg.realigned.bam")
execute(paste0(gatk.exec, " -T IndelRealigner -R ", reference.fa, " -targetIntervals ", indels.intervals , " -maxReads 150000 -I ", picard.rgfile, " -o ", gatk.realigned, " -known ", reference.indels, " -U ALLOW_N_CIGAR_READS"), gatk.realigned)

# Base Recalibration
gatk.cov1 <- paste0(filebase, ".1.covariates")
execute(paste0(gatk.exec, "-nct 4 -T BaseRecalibrator -R ", reference.fa.gz," -knownSites ", reference.snps, " -I ", gatk.realigned, " -o ", gatk.cov1," -U ALLOW_N_CIGAR_READS"), gatk.cov1)
gatk.recal <- paste0(filebase, ".aln.sort.rmdup.rg.realigned.recal.bam")
execute(paste0(gatk.exec, "-nct 4 -T PrintReads -R ", reference.fa.gz," -I ",gatk.realigned," -BQSR ", gatk.cov1," -U ALLOW_N_CIGAR_READS -o ", gatk.recal), gatk.recal)
gatk.cov2 <- paste0(filebase, ".2.covariates")
execute(paste0(gatk.exec, "-nct 4 -T BaseRecalibrator -R ", reference.fa.gz," -knownSites ", reference.snps, "  -I ",gatk.recal," -o ", gatk.cov2," -U ALLOW_N_CIGAR_READS"), gatk.cov2)
gatk.plots <- paste0(filebase, ".plots.pdf")
execute(paste0(gatk.exec, "-T AnalyzeCovariates -R ", reference.fa.gz," -before ", gatk.cov1, " -after ", gatk.cov2," -U ALLOW_N_CIGAR_READS -plots ", gatk.plots), gatk.plots)

# Done with the whole RNA-Seq pipeline
q("no")
