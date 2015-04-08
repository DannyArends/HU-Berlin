# Pipeline for DNA re-seq analysis on chicken
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

cmdlineargs   <- commandArgs(trailingOnly = TRUE)
fileBase      <- as.character(cmdlineargs[1])

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

referenceDir  <- "genomes"
reference     <- paste0(referenceDir, "/Leptin_Mouse.fasta")

## Do not forget to index the fasta file
#/opt/bwa-0.7.10/bwa index Gallus_gallus.Galgal4.fasta
## Also index the fasta for GATK
#java -jar /opt/picard-tools-1.99/CreateSequenceDictionary.jar R=Gallus_gallus.Galgal4.fasta O=Gallus_gallus.Galgal4.dict
#samtools faidx Gallus_gallus.Galgal4.fasta

########### ANALYSIS ###########

### Trimmomatic: Remove adapters and trim reads based on quality scores (1 to 2 hours) ###
logfile       <- paste0(fileBase,"log.txt")
trimmomatic   <- "/opt/Trimmomatic-0.32/"
gatk          <- "/opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"
picard        <- "/opt/picard-tools-1.99/"

inputfiles  <- c(paste0(fileBase, "_1.fq.gz"), paste0(fileBase, "_2.fq.gz"))
outputfiles <- c(paste0(fileBase, "_1.P_trimmed.fastq.gz"), paste0(fileBase, "_1.U_trimmed.fastq.gz"), 
                 paste0(fileBase, "_2.P_trimmed.fastq.gz"), paste0(fileBase, "_2.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.32.jar PE")
params  <- paste0("ILLUMINACLIP:", trimmomatic, "adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

command <- paste(cmdBase, inputfiles[1], inputfiles[2], outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)
#execute(command)

### BWA: Alignment against genome (2 to 8 hours), we pipe the output to the next step ###
# -v 2  : Verbose level 2 only errors and warnings
# -t 6  : Number of threads to use with BWA
command     <- paste0("/opt/bwa-0.7.10/bwa mem -v 2 -t 6 -A 3 -B 2 -U 4 -O 2 -E 0 -T 10 -a ", reference," ", outputfiles[1], " ", outputfiles[3], " | ")

### Convert SAM to BAM (1 hour), we pipe the output to the next step ###
# -Sb : Input sam, output bam
command     <- paste0(command, "samtools view -Sb - | ")

### Sort the BAM (1 hour) ###
# -@ : Number of CPU cores to use
# -m : Memory per CPU core
outputSBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.bam")
command      <- paste0(command, "samtools sort -@ 4 -m 2G -o - ", fileBase, " > ", outputSBAM)
execute(command)

### Add a read group ###
outputSRGBAM <- paste0(fileBase, "P_trimmed.aligned.sorted.rg.bam")
IDcode       <- substr(fileBase, 8, 12)
command      <- paste0("java -Xmx4g -jar ", picard, "AddOrReplaceReadGroups.jar INPUT=", outputSBAM, " OUTPUT=", outputSRGBAM, " CREATE_INDEX=false RGID=", IDcode, " RGLB=LIB", IDcode, " RGPL=Illumina RGPU=X RGSM=", IDcode)
execute(command)

### Move the file with the read group over the previous file ###
command <- paste0("mv ", outputSRGBAM, " ", outputSBAM)
execute(command)

### Index the BAM file (10 minutes) ###
outputSBAI   <- paste0(fileBase, "P_trimmed.aligned.sorted.bai")
command      <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
execute(command)

### Mark duplicates, using the Picard tools (~ 30 minutes) ###
outputSBAID    <- paste0(fileBase, "P_trimmed.aligned.sorted.dedup.bam")
outputMetrics  <- paste0(fileBase, ".metrics.txt")
command        <- paste0("java -jar ", picard, "MarkDuplicates.jar INPUT=", outputSBAM, " OUTPUT=", outputSBAID," METRICS_FILE=", outputMetrics," MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000")
execute(command)

### Index the BAM file (10 minutes) ###
outputSBAIDI <- paste0(fileBase, "P_trimmed.aligned.sorted.dedup.bai")
command      <- paste0("samtools index ", outputSBAID, " ", outputSBAIDI)
execute(command)

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSBAID)
execute(command)

### Indel Realign
outputSNPS    <- "output.snp.intervals"
outputSIBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.bam")
knownindels   <- paste0(referenceDir, "/Gallus_gallus.reordered.vcf")                               # Reference, download from: ftp://ftp.ensembl.org/pub/release-78/variation/vcf/gallus_gallus/
if(!file.exists(outputSNPS)){
  command <- paste0("java -Xmx4g -jar ", gatk, " -nt 4 -T RealignerTargetCreator -R ", reference, " -known ", knownindels, " -o ", outputSNPS, " -U ALLOW_N_CIGAR_READS")
  execute(command)                                                           # Call the GATK RealignerTargetCreator, only need to be done because the knownSNPlist does not change
}
command <- paste0("java -Xmx4g -jar ", gatk, " -T IndelRealigner -R ", reference, " -targetIntervals ", outputSNPS, " -I ", outputSBAID, " -o ",outputSIBAM," -known ",knownindels, " --consensusDeterminationModel KNOWNS_ONLY --fix_misencoded_quality_scores")
execute(command)                                                             # Call the GATK IndelRealigner

### Base Recalibration
knownsnps     <- paste0(referenceDir, "/Gallus_gallus.reordered.vcf")                               # Reference, download from: ftp://ftp.ensembl.org/pub/release-78/variation/vcf/gallus_gallus/
covfile1      <- paste0(fileBase, ".1.covariates")
covfile2      <- paste0(fileBase, ".2.covariates")
plotfile      <- paste0(fileBase, "recalibration.pdf")
outputSIRBAM  <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam")

command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIBAM," -o ",  covfile1, " -U ALLOW_N_CIGAR_READS")
execute(command)                                                            # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T PrintReads -R ", reference," -I ", outputSIBAM," -BQSR ", covfile1, " -U ALLOW_N_CIGAR_READS -o ", outputSIRBAM)
execute(command)                                                            # Call the GATK PrintReads
command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIRBAM," -o ", covfile2, " -U ALLOW_N_CIGAR_READS")
execute(command)                                                            # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -T AnalyzeCovariates -R ", reference, " -before ", covfile1, " -after ", covfile2, " -U ALLOW_N_CIGAR_READS -plots ", plotfile)
execute(command)                                                            # Call the GATK AnalyzeCovariates


### SNP calling
s1 <- "sample_6929P_trimmed.aligned.sorted.realigned.recalibrated.bam"
s2 <- "sample_6954P_trimmed.aligned.sorted.realigned.recalibrated.bam"
s3 <- "sample_8425-2P_trimmed.aligned.sorted.realigned.recalibrated.bam"
command <- paste0("java -Xmx4g -jar ", gatk," -T HaplotypeCaller -R ", reference," -I ", s1," -I ", s2, " -I ", s3, " --dbsnp ",knownsnps," -stand_call_conf 30 -stand_emit_conf 10 -L ", outputSNPS," -o output.raw.snps.indels.vcf")
execute(command)

#q("no")
