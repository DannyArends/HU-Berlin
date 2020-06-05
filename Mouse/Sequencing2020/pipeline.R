# Test with 150 mb file: nohup Rscript pipeline.R 19 B6-5 /home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35MNDSXY_L1_1.fq.gz /home/arends/NAS/Mouse/DNA/Sequencing/B6_2020/fastq/M_B6_5_FDSW202341663-1r_H35MNDSXY_L1_2.fq.gz &> nohup19-B6-5.out & 
cmdlineargs <- commandArgs(trailingOnly = TRUE)
seqID <- as.character(cmdlineargs[1])          # sequenceID
sampleID <- as.character(cmdlineargs[2])       # sampleID
R1 <- as.character(cmdlineargs[3])             # fastq R1 file
R2 <- as.character(cmdlineargs[4])             # fastq R2 file

outputfolder <- paste0("/home/arends/NAS/Mouse/DNA/Sequencing/Alignment2020/", seqID , "/")
referencefolder <- "/home/arends/NAS/Mouse/Reference_Genomes/GRCm38_68/"
reference <- paste0(referencefolder, "/GRCm38_68.fa")
logfile <- paste0(outputfolder, "log.txt")
outputBASE <- paste0(outputfolder, sampleID)

cat("--------------------------------------------------------------------\n")
cat("seqID:", seqID, "\n")
cat("sampleID:", sampleID, "\n")
cat("R1:", R1, "\n")
cat("R2:", R2, "\n")
cat("outputfolder:", outputfolder, "\n")
cat("outputBASE:", outputBASE, "\n")
cat("logfile:", logfile, "\n")
cat("reference:", reference, "\n")
cat("--------------------------------------------------------------------\n")

dir.create(outputfolder, showWarnings = FALSE, recursive = TRUE)
cat("", file = logfile)
cat("Created output folder\n", file = logfile, append = TRUE)

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

### Trimmomatic: Remove adapters and trim reads based on quality scores (1 to 2 hours) ###
trimmomatic <- "/home/arends/Github/Trimmomatic-0.39/"
adapters <- "/home/arends/Github/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

trimmomaticout <- paste0(outputfolder, "trimmed/")
dir.create(trimmomaticout, showWarnings = FALSE, recursive = TRUE)
cat("Created Trimmomatic output folder\n", file = logfile, append = TRUE)

outputfiles <- c(paste0(trimmomaticout, sampleID, "1.P_trimmed.fastq.gz"), paste0(trimmomaticout, sampleID, "1.U_trimmed.fastq.gz"), 
                 paste0(trimmomaticout, sampleID, "2.P_trimmed.fastq.gz"), paste0(trimmomaticout, sampleID, "2.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.39.jar PE")
params  <- paste0("ILLUMINACLIP:", adapters, ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

command <- paste(cmdBase, R1, R2, outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)

cat("Execute:", command, "\n", file = logfile, append = TRUE)
execute(command)
cat("-----------------------------------\n")

# Alignment using BWA
command <- paste0("~/Github/bwa-0.7.17/bwa mem -v 1 -t 4 -a ", reference," ", outputfiles[1], " ", outputfiles[3], " | ")

### Convert SAM to BAM (1 hour), we pipe the output to the next step ###
# -Sb : Input sam, output bam
command <- paste0(command, "samtools view -Sb - | ")

### Sort the BAM (1 hour) ###
# -@ : Number of CPU cores to use
# -m : Memory per CPU core
outputSBAM <- paste0(outputBASE, "P_trimmed.aligned.sorted.bam")
command <- paste0(command, "samtools sort -@ 2 - > ", outputSBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Add a read group ###
picard <- "/home/arends/Github/picard-2.22.3/picard.jar"
outputSRGBAM <- paste0(outputBASE, "P_trimmed.aligned.sorted.rg.bam")
command      <- paste0("java -Xmx4g -jar ", picard, " AddOrReplaceReadGroups INPUT=", outputSBAM, " OUTPUT=", outputSRGBAM, " CREATE_INDEX=false RGID=", seqID, " RGLB=LIB-", paste0(sampleID, "-",seqID), " RGPL=Illumina RGPU=X RGSM=", sampleID)
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

### Base Recalibration
gatk <- "/home/arends/Github/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar"

### Base Recalibration
knownsnps <- paste0(referencefolder, "/mgp.v5.merged.snps_all.dbSNP142.vcf.gz")
covfile <- paste0(outputBASE, ".covariates")
plotfile <- paste0(outputBASE, "recalibration.pdf")
outputSIRBAM <- paste0(outputBASE, "P_trimmed.aligned.sorted.dedup.recalibrated.bam")

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
outputSBAIDI <- paste0(outputBASE, "P_trimmed.aligned.sorted.dedup.recalibrated.bai")
command      <- paste0("samtools index ", outputSIRBAM, " ", outputSBAIDI)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSIRBAM)
cat("Execute:", command, "\n")
execute(command)
cat("-----------------------------------\n")

q("no")
