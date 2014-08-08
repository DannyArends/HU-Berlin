# Pipeline for RNA seq analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

cmdlineargs   <- commandArgs(trailingOnly = TRUE)
fileBase      <- as.character(cmdlineargs[1])

referenceDir  <- "genomes"
referenceName <- "Mus_musculus.GRCm38.74.dna"

########### SETUP THE REFERENCE AND KNOWN SNPS AND INDELS ###########

### Download the reference genome, and create a single big FASTA
reference     <- paste0(referenceDir, "/", referenceName, ".fasta")
#chromosomes <- c(1:19, "X", "Y", "MT")
#cat("", file=reference)
#for(chr in chromosomes){
#  system(paste0("wget -P temp ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa.gz"))
#  system(paste0("gunzip temp/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa.gz"))
#  fastadata <- readLines(paste0("temp/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa"))
#  cat(fastadata, sep="\n", file = reference, append = TRUE)
#  system(paste0("rm temp/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa.gz")) # Delete the temp fasta file
#}

### Download the known INDELS and SNPs
# system("wget -P genomes ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf.gz")
# system("wget -P genomes ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/mgp.v3.snps.rsIDdbSNPv137.vcf.gz")

########### PRE-PROCESSING THE REFERENCE ###########

### Create A BWA index file from the reference genome / transcriptome (10 to 100 minutes) ###
#system(paste0("/opt/bwa-0.7.10/bwa index ", reference))

### Create A SAMTOOLS index file from the reference genome ###
#system(paste0("samtools faidx ", reference))

### Create a Picard dictionary file
#referenceDict <- paste0(referenceDir, "/", referenceName, ".dict")
#system(paste0("java -jar /opt/picard-tools-1.99/CreateSequenceDictionary.jar R=", reference, " O=", referenceDict))

########### ANALYSIS ###########

### Trimmomatic: Remove adapters and trim reads based on quality scores (1 to 2 hours) ###
logfile       <- paste0(fileBase,"log.txt")
trimmomatic   <- "/opt/Trimmomatic-0.32/"
gatk          <- "/opt/GenomeAnalysisTK-2.7-2/GenomeAnalysisTK.jar"

cat("", file=logfile)     # Clear the logfile

inputfiles  <- c(paste0(fileBase, "R1_001.fastq.gz"), paste0(fileBase, "R2_001.fastq.gz"))
outputfiles <- c(paste0(fileBase, "R1_001.P_trimmed.fastq.gz"), paste0(fileBase, "R1_001.U_trimmed.fastq.gz"), 
                 paste0(fileBase, "R2_001.P_trimmed.fastq.gz"), paste0(fileBase, "R2_001.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.32.jar PE")
params  <- paste0("ILLUMINACLIP:", trimmomatic, "adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")

command <- paste(cmdBase, inputfiles[1], inputfiles[2], outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n")

### BWA: Alignment against genome (2 to 8 hours), we pipe the output to the next step ###
# -v 2  : Verbose level 2 only errors and warnings
# -t 6  : Number of threads to use with BWA
command     <- paste0("/opt/bwa-0.7.10/bwa mem -v 2 -t 6 ", reference," ", outputfiles[1], " ", outputfiles[3], " | ")

### Convert SAM to BAM (1 hour), we pipe the output to the next step ###
# -Sb : Input sam, output bam
command     <- paste0(command, "samtools view -Sb - | ")

### Sort the BAM (1 hour) ###
# -@ : Number of CPU cores to use
# -m : Memory per CPU core
filePrefix   <- fileBase
outputSBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.bam")
command      <- paste0(command, "samtools sort -@ 4 -m 2G -o - ", filePrefix, " > ", outputSBAM)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))

### Add a read group ###
outputSRGBAM <- paste0(fileBase, "P_trimmed.aligned.sorted.rg.bam")
IDcode       <- substr(fileBase, 1, 4)
command      <- paste0("java -Xmx4g -jar /opt/picard-tools-1.99/AddOrReplaceReadGroups.jar INPUT=", outputSBAM, " OUTPUT=", outputSRGBAM," CREATE_INDEX=false RGID=",IDcode," RGLB=LIB",IDcode," RGPL=Illumina RGPU=X RGSM=", IDcode)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))

### Move the file with the read group over the previous file
command <- paste0("mv ", outputSRGBAM, " ", outputSBAM)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))

### Index the BAM file (10 minutes) ###
outputSBAI   <- paste0(fileBase, "P_trimmed.aligned.sorted.bai")
command      <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSBAM)
cat(system(command, intern=TRUE), file=paste0(fileBase,"stats.txt"), sep="\n"))

### Indel Realign
outputSNPS    <- "output.snp.intervals"
outputSIBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.bam")
knownindels   <- "genomes/mgp.v3.indels.rsIDdbSNPv137.vcf"                              # Reference, download from: ftp://ftp-mouse.sanger.ac.uk/
if(!file.exists(outputSNPS)){
  command <- paste0("java -Xmx4g -jar ", gatk, " -nt 4 -T RealignerTargetCreator -R ", reference, " -known ", knownindels, " -o ", outputSNPS, " -U ALLOW_N_CIGAR_READS")
  cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))   # Call the GATK RealignerTargetCreator, only need to be done because the knownSNPlist does not change
}
command <- paste0("java -Xmx4g -jar ", gatk, " -T IndelRealigner -R ", reference, " -targetIntervals ", outputSNPS, " -maxReads 150000 -I ", outputSBAM, " -o ",outputSIBAM," -known ",knownindels, " -U ALLOW_N_CIGAR_READS")
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))     # Call the GATK IndelRealigner

### Base Recalibration
knownsnps     <- "genomes/mgp.v3.snps.rsIDdbSNPv137.vcf"                                # Reference, download from: ftp://ftp-mouse.sanger.ac.uk/
covfile1      <- paste0(fileBase, "1.covariates")
covfile2      <- paste0(fileBase, "2.covariates")
plotfile      <- paste0(fileBase, "recalibration.pdf")
outputSIRBAM  <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam")

command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIBAM," -o ",  covfile1, " -U ALLOW_N_CIGAR_READS")
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))   # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T PrintReads -R ", reference," -I ", outputSIBAM," -BQSR ", covfile1, " -U ALLOW_N_CIGAR_READS -o ", outputSIRBAM)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))   # Call the GATK PrintReads
command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIRBAM," -o ", covfile2, " -U ALLOW_N_CIGAR_READS")
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))   # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -T AnalyzeCovariates -R ", reference, " -before ", covfile1, " -after ", covfile2, " -U ALLOW_N_CIGAR_READS -plots ", plotfile)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))   # Call the GATK AnalyzeCovariates

### SNP calling -> VCF for NoverlSNPer
# See http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_haplotypecaller_HaplotypeCaller.html
outputVCF    <- paste0(fileBase, ".snps.vcf")
settings     <- "-recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0" 

command <- paste0("java -Xmx4g -jar ", gatk, " -T HaplotypeCaller -R ", reference, " -I ", outputSIRBAM, "  --dbsnp ", knownsnps, " ", settings, " -o ", outputVCF)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n"))   # Call the GATK HaplotypeCaller

### Expression via Express, DeSeq2
#library("GenomicAlignments")
#library("Rsamtools")

#hse         <- makeTranscriptDbFromGFF("/path/to/your/genemodel.GTF", format = "gtf")
#exonsByGene <- exonsBy(hse, by = "gene")
#files       <- list.files("/path/to/bam/files", pattern="bam$", full=TRUE)
#bamfiles    <- BamFileList(files, yieldSize=100000)
#se          <- summarizeOverlaps(exonsByGene, bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

#head(assay(se))
q("no")