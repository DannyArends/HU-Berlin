# Pipeline for DNA re-seq analysis on mouse
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Mar, 2019
# first written Jan, 2015

#
# zcat /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ext_L7258-1_BFMI860-S12_S1_R1_001.fastq.gz | head -n 400000 > /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/subset_R1_001.fastq
# zcat /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ext_L7258-1_BFMI860-S12_S1_R2_001.fastq.gz | head -n 400000 > /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/subset_R2_001.fastq
#

# nohup Rscript pipeline_Mar19.R subset_R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ /halde/BFMI_Alignment_Mar19/ BFMI860-12 &> subset_nohup.out &

# nohup Rscript pipeline_Mar19.R ext_L7258-1_BFMI860-S12_S1_R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ /halde/BFMI_Alignment_Mar19/ BFMI860-12 &> ext_L7258-1_BFMI860-S12_S1_nohup.out &
# nohup Rscript pipeline_Mar19.R ext_L7258-1_BFMI860-S12_S9_R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ /halde/BFMI_Alignment_Mar19/ BFMI860-12 &> ext_L7258-1_BFMI860-S12_S9_nohup.out &
# nohup Rscript pipeline_Mar19.R ext_L7258-2_BFMI860-S12_S2_R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ /halde/BFMI_Alignment_Mar19/ BFMI860-12 &> ext_L7258-2_BFMI860-S12_S2_nohup.out &
# nohup Rscript pipeline_Mar19.R ext_L7258-1_BFMI860-S12_S23_R /home/danny/NAS/Mouse/DNA/Sequencing/BFMI860_Feb2019/ /halde/BFMI_Alignment_Mar19/ BFMI860-12 &> ext_L7258-1_BFMI860-S12_S23_nohup.out &

# Afterwards merge them into a single BAM: samtools merge -@ 8 merged.bam ext_L7258-1_BFMI860-S12_S1_RP_trimmed.aligned.sorted.dedup.bam ext_L7258-1_BFMI860-S12_S23_RP_trimmed.aligned.sorted.dedup.bam ext_L7258-1_BFMI860-S12_S9_RP_trimmed.aligned.sorted.dedup.bam ext_L7258-2_BFMI860-S12_S2_RP_trimmed.aligned.sorted.dedup.bam
# Sort the merged BAM: samtools sort -@ 8 merged.bam > merged_sorted.bam
# Index the merged sorted BAM: samtools index -@ 8 merged_sorted.bam merged_sorted.bai


cmdlineargs <- commandArgs(trailingOnly = TRUE)
inBaseName <- as.character(cmdlineargs[1])              # File base of the fastq files
inputFolder <- as.character(cmdlineargs[2])             # Input folder where the fastq files are found
outFolder <- as.character(cmdlineargs[3])               # Output folder for the script
IDcode <- as.character(cmdlineargs[4])                  # Read group ID

inputBASE <- paste0(inputFolder, inBaseName)
outputBASE <- paste0(outFolder, inBaseName)

cat("--------------------------------------------------------------------\n")
cat("inBaseName:", inBaseName, "\n")
cat("inputFolder:", inputFolder, "\n")
cat("inputBASE:", inputBASE, "\n")
cat("outFolder:", outFolder, "\n")
cat("outputBASE:", outputBASE, "\n")
cat("IDcode:", IDcode, "\n")
cat("--------------------------------------------------------------------\n")

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

referenceDir  <- "/home/danny/References/Mouse/GRCm38_95"
reference     <- paste0(referenceDir, "/Mus_musculus.GRCm38.dna.toplevel.fa.gz")

### Trimmomatic: Remove adapters and trim reads based on quality scores (1 to 2 hours) ###
logfile       <- paste0(outputBASE, "log.txt")
trimmomatic   <- "/home/danny/Github/trimmomatic/trimmomatic-0.38/dist/jar/"
adapters      <- "/home/danny/Github/trimmomatic/trimmomatic-0.38/adapters/TruSeq3-PE.fa"
gatk          <- "/home/danny/Github/gatk-4.1.0.0/gatk-package-4.1.0.0-local.jar"
picard        <- "/home/danny/Github/picard-2.19.0/picard.jar"

inputfiles  <- c(paste0(inputBASE, "1_001.fastq.gz"), paste0(inputBASE, "2_001.fastq.gz"))
outputfiles <- c(paste0(outputBASE, "1.P_trimmed.fastq.gz"), paste0(outputBASE, "1.U_trimmed.fastq.gz"), 
                 paste0(outputBASE, "2.P_trimmed.fastq.gz"), paste0(outputBASE, "2.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.38.jar PE")
params  <- paste0("ILLUMINACLIP:", adapters, ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

command <- paste(cmdBase, inputfiles[1], inputfiles[2], outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)
cat("-------------------------------------------------------------------- Trimmomatic Start\n")
execute(command)
cat("-------------------------------------------------------------------- Trimmomatic Done\n")

# Alignment using BWA
command     <- paste0("~/Github/bwa/bwa mem -v 2 -t 12 -a ", reference," ", outputfiles[1], " ", outputfiles[3], " | ")

### Convert SAM to BAM (1 hour), we pipe the output to the next step ###
# -Sb : Input sam, output bam
command     <- paste0(command, "samtools view -Sb - | ")

### Sort the BAM (1 hour) ###
# -@ : Number of CPU cores to use
# -m : Memory per CPU core
outputSBAM   <- paste0(outputBASE, "P_trimmed.aligned.sorted.bam")
command      <- paste0(command, "samtools sort -@ 8 -m 4G - > ", outputSBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Add a read group ###
outputSRGBAM <- paste0(outputBASE, "P_trimmed.aligned.sorted.rg.bam")
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
outputSBAI   <- paste0(outputBASE, "P_trimmed.aligned.sorted.bai")
command      <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Mark duplicates, using the Picard tools (~ 30 minutes) ###
outputSBAID    <- paste0(outputBASE, "P_trimmed.aligned.sorted.dedup.bam")
outputMetrics  <- paste0(outputBASE, ".metrics.txt")
command        <- paste0("java -jar ", picard, " MarkDuplicates INPUT=", outputSBAM, " OUTPUT=", outputSBAID," METRICS_FILE=", outputMetrics," MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000")
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Index the BAM file (10 minutes) ###
outputSBAIDI <- paste0(outputBASE, "P_trimmed.aligned.sorted.dedup.bai")
command      <- paste0("samtools index ", outputSBAID, " ", outputSBAIDI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSBAID)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

q("no")
