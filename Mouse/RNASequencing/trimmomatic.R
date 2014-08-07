# Perform read trimming using Trimmomatic v0.32
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

### Download the known snps, and convert to SNPs
# wget ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf.gz
# mgp.v3.indels.rsIDdbSNPv137.vcf.gz

########### PRE-PROCESSING THE REFERENCE ###########

### Create A BWA index file from the reference genome / transcriptome (10 to 100 minutes) ###
# Note this has to be performed only once, the index file can be reused

#command <- "/opt/bwa-0.7.10/bwa index Mus_musculus.GRCm38.74.dna.fasta"
#system(command)

### Create A SAMTOOLS index file from the reference genome ###
#command <- "samtools faidx Mus_musculus.GRCm38.74.dna.fasta"
#system(command)

### Create a picard dictionary file
#command <- paste0("java -jar /opt/picard-tools-1.99/CreateSequenceDictionary.jar R=Mus_musculus.GRCm38.74.dna.fasta O=Mus_musculus.GRCm38.74.dna.dict")


########### ANALYSIS ###########

### Trimmomatic: Remove adapters and trim reads based on quality scores (1 to 2 hours) ###
reference     <- "genomes/Mus_musculus.GRCm38.74.dna.fasta"
fileBase      <- "4868_GCCAAT_L001_"
trimmomatic   <- "/opt/Trimmomatic-0.32/"
gatk          <- "/opt/GenomeAnalysisTK-2.7-2/GenomeAnalysisTK.jar"

inputfiles  <- c(paste0(fileBase, "R1_001_S.fastq.gz"), paste0(fileBase, "R2_001_S.fastq.gz"))
outputfiles <- c(paste0(fileBase, "R1_001_S.P_trimmed.fastq.gz"), paste0(fileBase, "R1_001_S.U_trimmed.fastq.gz"), 
                 paste0(fileBase, "R2_001_S.P_trimmed.fastq.gz"), paste0(fileBase, "R2_001_S.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.32.jar PE")
params  <- paste0("ILLUMINACLIP:", trimmomatic, "adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

command <- paste(cmdBase, inputfiles[1], inputfiles[2], outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)
system(command)

### BWA: Alignment against genome (2 to 8 hours), we pipe the output to the next step ###
# -v 2  : Verbose level 2 only errors and warnings
# -t 6  : Number of threads to use with BWA
command     <- paste0("/opt/bwa-0.7.10/bwa mem -v 2 -t 6 ", reference," ", outputfiles[1], " ", outputfiles[3], " | ")

### Convert SAM to BAM (1 hour), we pipe the output to the next step ###
# -Sb : Input sam, output bam
command     <- paste0(command, "samtools view -Sb - | ")

### Sort the BAM (1 hour) ###
filePrefix   <- fileBase
outputSBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.bam")
command      <- paste0(command, "samtools sort -@ 4 -m 2G -o - ", filePrefix, " > ", outputSBAM)
system(command)

### Add a read group ###
outputSRGBAM <- paste0(fileBase, "P_trimmed.aligned.sorted.rg.bam")
IDcode       <- substr(fileBase, 1, 4)
command      <- paste0("java -Xmx4g -jar /opt/picard-tools-1.99/AddOrReplaceReadGroups.jar INPUT=", outputSBAM, " OUTPUT=", outputSRGBAM," CREATE_INDEX=true RGID=",IDcode," RGLB=LIB",IDcode," RGPL=Illumina RGPU=X RGSM=", IDcode)
system(command)

### Move the file with the read group over the previous file
command <- paste0("mv ", outputSRGBAM, " ", outputSBAM)
system(command)

### Index the BAM file (10 minutes) ###
outputSBAI   <- paste0(fileBase, "P_trimmed.aligned.sorted.bai")
command      <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
system(command)

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSBAM)
system(command)

### Indel Realign
outputSNPS    <- paste0(fileBase, ".snps.list")
outputSIBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.bam")

command <- paste0("java -Xmx4g -jar ", gatk, " -T RealignerTargetCreator -R ", reference, " -I ", outputSBAM, " -o ", outputSNPS, " -U ALLOW_N_CIGAR_READS")
system(command)   # Call the GATK RealignerTargetCreator
command <- paste0("java -Xmx4g -jar ", gatk, " -T IndelRealigner -R ", reference, " -targetIntervals ", outputSNPS, " -maxReads 150000 -I ", outputSBAM, " -o ",outputSIBAM," -U ALLOW_N_CIGAR_READS")
system(command)   # Call the GATK IndelRealigner

### Expression via Express, DeSeq2
#library("GenomicAlignments")
#library("Rsamtools")

#hse         <- makeTranscriptDbFromGFF("/path/to/your/genemodel.GTF", format = "gtf")
#exonsByGene <- exonsBy(hse, by = "gene")
#files       <- list.files("/path/to/bam/files", pattern="bam$", full=TRUE)
#bamfiles    <- BamFileList(files, yieldSize=100000)
#se          <- summarizeOverlaps(exonsByGene, bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

#head(assay(se))

### Base Recalibration
snpome        <- "genomes/mgp.v3.indels.rsIDdbSNPv137.vcf"
covfile1      <- paste0(fileBase, ".1.covariates")
covfile2      <- paste0(fileBase, ".2.covariates")
plotfile      <- paste0(fileBase, ".recalibration.pdf")
outputSIRBAM  <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam")

command <- paste0("java -Xmx4g -jar ", gatk, " -T BaseRecalibrator -R ", reference, " -knownSites ", snpome, " -I ", outputSIBAM," -o ",  covfile1, " -U ALLOW_N_CIGAR_READS")
system(command)   # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -T PrintReads -R ", reference," -I ", outputSIBAM," -BQSR ", covfile, " -U ALLOW_N_CIGAR_READS -o ", outputSIRBAM)
system(command)   # Call the GATK PrintReads
command <- paste0("java -Xmx4g -jar ", gatk, " -T BaseRecalibrator -R ", reference, " -knownSites ", snpome, " -I ", outputSIRBAM," -o ", covfile2, " -U ALLOW_N_CIGAR_READS")
system(command)   # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -T AnalyzeCovariates -R ", reference, " -before ", covfile1, " -after ", covfile2, " -U ALLOW_N_CIGAR_READS -plots ", plotfile)
system(command)   # Call the GATK AnalyzeCovariates

### Haplotype calling -> VCF for NoverlSNPer
outputVCF    <- paste0(fileBase, ".snps.vcf")
settings     <- "-stand_call_conf 30.0 -stand_emit_conf 10.0"

command <- paste0("java -Xmx4g -jar ", gatk, " -T HaplotypeCaller -R ", reference, " -I ", outputSIRBAM, "  --dbsnp ", snpome," ", settings," -o ", outputVCF, " -U ALLOW_N_CIGAR_READS")
system(command)   # Call the GATK HaplotypeCaller

