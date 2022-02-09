#
# zcat /home/danny/NAS/Mouse/RNAseqKourosh/AKR-1739/HF3W7DRXY/L001/AKR-1739_S1_L001_R1_001.fastq.gz | head -n 100
# zcat /home/danny/NAS/Mouse/RNAseqKourosh/alignments/trimmed/AKR-1739trimmed.fastq.gz | head -n 100
#

fastq <- c(
  "AKR-1739/HF3W7DRXY/L001/AKR-1739_S1_L001_R1_001.fastq.gz",
  "AKR-1749/HF3W7DRXY/L001/AKR-1749_S2_L001_R1_001.fastq.gz",
  "AKR-1756/HF3W7DRXY/L001/AKR-1756_S3_L001_R1_001.fastq.gz",
  "AKR-1761/HF3W7DRXY/L001/AKR-1761_S4_L001_R1_001.fastq.gz",
  
  "B6-1_1009437/HF3W7DRXY/L001/B6-1_1009437_S5_L001_R1_001.fastq.gz",
  "B6-2_1009438/HF3W7DRXY/L001/B6-2_1009438_S6_L001_R1_001.fastq.gz",
  "B6-3_1009442/HF3W7DRXY/L001/B6-3_1009442_S7_L001_R1_001.fastq.gz",
  "B6-4_1009443/HF3W7DRXY/L001/B6-4_1009443_S8_L001_R1_001.fastq.gz",
  
  "BFMI-1_860-29698/HF3W7DRXY/L001/BFMI-1_860-29698_S9_L001_R1_001.fastq.gz",
  "BFMI-2_860-29699/HF3W7DRXY/L001/BFMI-2_860-29699_S10_L001_R1_001.fastq.gz",
  "BFMI-3_860-29700/HF3W7DRXY/L001/BFMI-3_860-29700_S11_L001_R1_001.fastq.gz",
  "BFMI-4_860-29692/HF3W7DRXY/L001/BFMI-4_860-29692_S12_L001_R1_001.fastq.gz",
  
  "SJL-1216/HF3W7DRXY/L001/SJL-1216_S13_L001_R1_001.fastq.gz",
  "SJL-1223/HF3W7DRXY/L001/SJL-1223_S14_L001_R1_001.fastq.gz",
  "SJL-1227/HF3W7DRXY/L001/SJL-1227_S15_L001_R1_001.fastq.gz",
  "SJL-1234/HF3W7DRXY/L001/SJL-1234_S16_L001_R1_001.fastq.gz"
)
for(sampleID in 2:length(fastq)){
fileloc <- "/home/danny/NAS/Mouse/RNAseqKourosh/"
sampleName <- unlist(lapply(strsplit(fastq[sampleID], "/"), "[",1))
R1 <- paste0(fileloc, fastq[sampleID])

outputfolder <- "/home/danny/NAS/Mouse/RNAseqKourosh/alignments/"
logfile <- paste0(outputfolder, "log", sampleName,".txt")
referencefolder <- "/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68"
reference <- paste0(referencefolder, "/GRCm38_68.fa")

outputBASE <- paste0(outputfolder, sampleName)

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

### cutadapt: Remove adapters and trim reads ###

trimmomaticout <- paste0(outputfolder, "trimmed/")
dir.create(trimmomaticout, showWarnings = FALSE, recursive = TRUE)
cat("Created Trimmomatic output folder\n", file = logfile, append = TRUE)

outputfiles <- c(paste0(trimmomaticout, sampleName, ".t1.fastq.gz"), paste0(trimmomaticout, sampleName, ".t2.fastq.gz"))

command <- paste0("cutadapt -a TGGAATTCTCGGGTGCCAAGG --minimum-length 23 ", R1, " -o ", outputfiles[1])
cat("Execute:", command, "\n")
#execute(command)
command <- paste0("cutadapt -u 4 -u -4 ", R1, " -o ", outputfiles[2])
cat("Execute:", command, "\n")
#execute(command)

cat("-----------------------------------\n")


# Alignment using BWA
command <- paste0("/home/danny/Github/STAR/source/STAR --runMode alignReads --readFilesCommand zcat --genomeDir=/home/danny/NAS/Mouse/Reference_Genomes/GRCm38_68/STAR --readFilesIn ", outputfiles[2], " --outSAMtype SAM --outStd SAM | ")

### Convert SAM to BAM (1 hour), we pipe the output to the next step ###
# -Sb : Input sam, output bam
command <- paste0(command, "samtools view -Sb - | ")

### Sort the BAM (1 hour) ###
# -@ : Number of CPU cores to use
# -m : Memory per CPU core
outputSBAM <- paste0(outputBASE, ".STAR.trimmed.aligned.sorted.bam")
command <- paste0(command, "samtools sort -@ 2 - > ", outputSBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Add a read group ###
picard <- "/home/danny/Github/picard-2.25.6/picard.jar"
outputSRGBAM <- paste0(outputBASE, ".STAR.trimmed.aligned.sorted.rg.bam")
command      <- paste0("java -Xmx4g -jar ", picard, " AddOrReplaceReadGroups INPUT=", outputSBAM, " OUTPUT=", outputSRGBAM, " CREATE_INDEX=false RGID=", sampleID, " RGLB=LIB-", sampleName, " RGPL=Illumina RGPU=X RGSM=", sampleID)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Move the file with the read group over the previous file ###
command <- paste0("mv ", outputSRGBAM, " ", outputSBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Index the BAM file (10 minutes) ###
outputSBAI   <- paste0(outputBASE, ".STAR.trimmed.aligned.sorted.bai")
command      <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Mark duplicates, using the Picard tools (~ 30 minutes) ###
outputSBAID    <- paste0(outputBASE, ".STAR.trimmed.aligned.sorted.dedup.bam")
outputMetrics  <- paste0(outputBASE, ".STAR.metrics.txt")
command        <- paste0("java -jar ", picard, " MarkDuplicates INPUT=", outputSBAM, " OUTPUT=", outputSBAID," METRICS_FILE=", outputMetrics," MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000")
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Index the BAM file (10 minutes) ###
outputSBAIDI <- paste0(outputBASE, ".STAR.trimmed.aligned.sorted.dedup.bai")
command      <- paste0("samtools index ", outputSBAID, " ", outputSBAIDI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Base Recalibration
gatk <- "/home/danny/Github/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar"

### Base Recalibration
knownsnps <- paste0(referencefolder, "/mgp.v5.merged.snps_all.dbSNP142.vcf.gz")
covfile <- paste0(outputBASE, ".STAR.covariates")
plotfile <- paste0(outputBASE, ".STAR.recalibration.pdf")
outputSIRBAM <- paste0(outputBASE, ".STAR.trimmed.aligned.sorted.dedup.recalibrated.bam")

command <- paste0("java -Xmx4g -jar ", gatk, " BaseRecalibrator -R ", reference, " --known-sites ", knownsnps, " -I ", outputSBAID," -O ",  covfile)
cat("Execute:", command, "\n")  # Call the GATK BaseRecalibrator
execute(command)
cat("-----------------------------------\n")

command <- paste0("java -Xmx4g -jar ", gatk, " ApplyBQSR -R ", reference," -I ", outputSBAID," --bqsr-recal-file ", covfile, " -O ", outputSIRBAM)
cat("Execute:", command, "\n")  # Call the GATK ApplyBQSR
execute(command)
cat("-----------------------------------\n")

command <- paste0("java -Xmx4g -jar ", gatk, " AnalyzeCovariates -bqsr ", covfile, " -plots ", plotfile)
cat("Execute:", command, "\n")  # Call the GATK AnalyzeCovariates
execute(command)
cat("-----------------------------------\n")

### Index the BAM file (10 minutes) ###
outputSBAIDI <- paste0(outputBASE, ".STAR.trimmed.aligned.sorted.dedup.recalibrated.bai")
command      <- paste0("samtools index ", outputSIRBAM, " ", outputSBAIDI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSIRBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")
}
q("no")
