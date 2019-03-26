# Pipeline for DNA re-seq analysis on mouse
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Mar, 2019
# first written Jan, 2015

cmdlineargs   <- commandArgs(trailingOnly = TRUE)
inputBASE     <- as.character(cmdlineargs[1])                        # File base of the fastq files
inName        <- unlist(strsplit(inputBASE, "/"))
inName        <- gsub(".bam", "", inName[length(inName)])

outFolder     <- as.character(cmdlineargs[2])                         # Output folder for the script
fileBase      <- paste0(outFolder, "/", inName)

IDcode       <- as.character(cmdlineargs[3])                          # Read group ID

cat("inputBASE:", inputBASE, "\n")
cat("inName:", inName, "\n")
cat("fileBase:", fileBase, "\n")

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

referenceDir  <- "~/References/Mouse/GRCm38_95"
reference     <- paste0(referenceDir, "/Mus_musculus.GRCm38.dna.toplevel.fa.gz")

### Trimmomatic: Remove adapters and trim reads based on quality scores (1 to 2 hours) ###
logfile       <- paste0(fileBase, "log.txt")
trimmomatic   <- "/home/danny/Github/trimmomatic/trimmomatic-0.38/dist/jar/"
adapters      <- "/home/danny/Github/trimmomatic/trimmomatic-0.38/adapters/TruSeq3-PE.fa"
gatk          <- "/home/danny/Github/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
picard        <- "/home/danny/Github/picard-2.19.0/picard.jar"

inputfiles  <- c(paste0(inputBASE, "1_001.fastq.gz"), paste0(inputBASE, "2_001.fastq.gz"))
outputfiles <- c(paste0(fileBase, "1.P_trimmed.fastq.gz"), paste0(fileBase, "1.U_trimmed.fastq.gz"), 
                 paste0(fileBase, "2.P_trimmed.fastq.gz"), paste0(fileBase, "2.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.38.jar PE")
params  <- paste0("ILLUMINACLIP:", adapters, ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

command <- paste(cmdBase, inputfiles[1], inputfiles[2], outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)
cat("Skip:", command, "\n")
#execute(command)

cat("-----------------------------------\n")

# Alignment using BWA
command     <- paste0("~/Github/bwa/bwa mem -v 2 -t 6 -A 3 -B 2 -U 4 -O 2 -E 0 -T 10 -a ", reference," ", outputfiles[1], " ", outputfiles[3], " | ")

### Convert SAM to BAM (1 hour), we pipe the output to the next step ###
# -Sb : Input sam, output bam
command     <- paste0(command, "samtools view -Sb - | ")

### Sort the BAM (1 hour) ###
# -@ : Number of CPU cores to use
# -m : Memory per CPU core
outputSBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.bam")
command      <- paste0(command, "samtools sort -@ 8 -m 4G - > ", outputSBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")
q("no")

### Add a read group ###
outputSRGBAM <- paste0(fileBase, "P_trimmed.aligned.sorted.rg.bam")
command      <- paste0("java -Xmx4g -jar ", picard, " AddOrReplaceReadGroups INPUT=", outputSBAM, " OUTPUT=", outputSRGBAM, " CREATE_INDEX=false RGID=", IDcode, " RGLB=LIB", IDcode, " RGPL=Illumina RGPU=X RGSM=", IDcode)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Move the file with the read group over the previous file ###
command <- paste0("mv ", outputSRGBAM, " ", outputSBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Index the BAM file (10 minutes) ###
outputSBAI   <- paste0(fileBase, "P_trimmed.aligned.sorted.bai")
command      <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Mark duplicates, using the Picard tools (~ 30 minutes) ###
outputSBAID    <- paste0(fileBase, "P_trimmed.aligned.sorted.dedup.bam")
outputMetrics  <- paste0(fileBase, ".metrics.txt")
command        <- paste0("java -jar ", picard, " MarkDuplicates INPUT=", outputSBAM, " OUTPUT=", outputSBAID," METRICS_FILE=", outputMetrics," MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000")
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Index the BAM file (10 minutes) ###
outputSBAIDI <- paste0(fileBase, "P_trimmed.aligned.sorted.dedup.bai")
command      <- paste0("samtools index ", outputSBAID, " ", outputSBAIDI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSBAID)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Indel Realign
outputSNPS    <- "output.snp.intervals"
outputSIBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.bam")
knownindels   <- paste0(referenceDir, "/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz")
if(!file.exists(outputSNPS)){
  command <- paste0("java -Xmx4g -jar ", gatk, " -nt 4 -T RealignerTargetCreator -R ", reference, " -known ", knownindels, " -o ", outputSNPS, " -U ALLOW_N_CIGAR_READS")
  cat("Execute:", command, "\n")
  execute(command)
  cat("-----------------------------------\n")
}
command <- paste0("java -Xmx4g -jar ", gatk, " -T IndelRealigner --filter_bases_not_stored -R ", reference, " -targetIntervals ", outputSNPS, " -I ", outputSBAID, " -o ",outputSIBAM," -known ",knownindels, " --consensusDeterminationModel KNOWNS_ONLY")
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")


### Base Recalibration
knownsnps     <- paste0(referenceDir, "/mgp.v5.merged.snps_all.dbSNP142.vcf.gz")
covfile1      <- paste0(fileBase, ".1.covariates")
covfile2      <- paste0(fileBase, ".2.covariates")
plotfile      <- paste0(fileBase, "recalibration.pdf")
outputSIRBAM  <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam")

command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIBAM," -o ",  covfile1, " -U ALLOW_N_CIGAR_READS")
cat("Execute:", command, "\n")  # Call the GATK BaseRecalibrator
execute(command)
cat("-----------------------------------\n")

command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T PrintReads -R ", reference," -I ", outputSIBAM," -BQSR ", covfile1, " -U ALLOW_N_CIGAR_READS -o ", outputSIRBAM)
cat("Execute:", command, "\n")  # Call the GATK PrintReads
execute(command)
cat("-----------------------------------\n")

command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIRBAM," -o ", covfile2, " -U ALLOW_N_CIGAR_READS")
cat("Execute:", command, "\n")  # Call the GATK BaseRecalibrator
execute(command)
cat("-----------------------------------\n")

command <- paste0("java -Xmx4g -jar ", gatk, " -T AnalyzeCovariates -R ", reference, " -before ", covfile1, " -after ", covfile2, " -U ALLOW_N_CIGAR_READS -plots ", plotfile)
cat("Execute:", command, "\n")  # Call the GATK AnalyzeCovariates
execute(command)
cat("-----------------------------------\n")

q("no")
