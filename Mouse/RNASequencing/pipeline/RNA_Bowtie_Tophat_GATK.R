## RNA Seq analysis pipeline, using bowtie2 and tophat2 for alignment, GATK for post-processing
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jun, 2015
# first written Jun, 2015

### Locations of tools
trimmomatic   <- "/opt/Trimmomatic-0.32/"
gatk          <- "/opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"

### Locations of directories
static.dir      <- "Reference/"; dir.create(static.dir)
input.dir       <- "FastQ/";  dir.create(input.dir)
output.dir      <- "Analysis/";  dir.create(output.dir)
log.dir         <- "Log/";  dir.create(log.dir)

### Reference files
reference.stem  <- "Mus_musculus.GRCm38"
reference.name  <- paste0(reference.stem, ".dna.toplevel.fa.gz")
reference.loc   <- paste0(static.dir, reference.name)
reference.unzip <- paste0(static.dir, paste0(reference.stem, ".dna.toplevel.fa"))
genesets.name   <- paste0(reference.stem, ".80.gtf.gz")
genesets.loc    <- paste0(static.dir, genesets.name)

### Read the sample description file
samples <- read.table("Experiment/SampleDescription.txt", sep="\t", header=TRUE)
sample.names <- apply(samples[, c("Lib_id","TagLane")], 1, paste0, collapse="_")
names(sample.names) <- samples[, "Lib_id"]

### Download the files we need (reference, genesets)
if(!file.exists(reference.loc)) download.file(paste0("ftp://ftp.ensembl.org/pub/release-80/fasta/mus_musculus/dna/", reference.name), reference.loc)
if(!file.exists(genesets.loc))  download.file(paste0("ftp://ftp.ensembl.org/pub/release-80/gtf/mus_musculus/", genesets.name), genesets.loc)

### Unzip the reference genome for bowtie2
command <- paste0("gunzip -c ", reference.loc, " > ", reference.unzip)
if(!file.exists(reference.unzip)) system(command, intern = TRUE)

### Build the bowtie2 index file for the reference genome
bowtie.index <- paste0(static.dir, reference.stem, ".dna.bt2idx")

command <- paste0("bowtie2-build -f ", reference.unzip," ", bowtie.index)
if(!file.exists(paste0(bowtie.index,".1.bt2"))) system(command, intern = TRUE)

### TODO: Get rid of the bed READs

### Trim RNA-Seq read files using Trimmomatic
sample.id    <- as.character("4423")
sample.name  <- sample.names[sample.id]
inputfiles   <- c(paste0(input.dir, sample.name, "_R1_001.subset.fastq.gz"), paste0(input.dir, sample.name, "_R2_001.subset.fastq.gz"))
outputfiles  <- c(paste0(output.dir, sample.name, "_R1_001.P_trimmed.fastq.gz"), paste0(output.dir, sample.name, "_R1_001.U_trimmed.fastq.gz"), 
                  paste0(output.dir, sample.name, "_R2_001.P_trimmed.fastq.gz"), paste0(output.dir, sample.name, "_R2_001.U_trimmed.fastq.gz"))

cmdBase <- paste0("java -jar ", trimmomatic, "trimmomatic-0.32.jar PE")
params  <- paste0("ILLUMINACLIP:", trimmomatic, "adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36")
command <- paste(cmdBase, inputfiles[1], inputfiles[2], outputfiles[1], outputfiles[2], outputfiles[3], outputfiles[4], params)
if(!file.exists(outputfiles[1])) system(command, intern = TRUE)

### Execute and run Tophat2 to align the trimmed files
sample.loc      <- paste0(input.dir,  sample.name)
log.file        <- paste0(log.dir,    sample.name, ".th2.log")
output.aligned  <- paste0(output.dir, sample.name, "_trim_align.bam")

th2.rg       <- paste0("--rg-id ", sample.id," --rg-sample ", sample.name, " --rg-library RNA-seq --rg-platform Illumina")
th2.options  <- paste0(th2.rg, " --prefilter-multihits --b2-sensitive --read-gap-length 12 --read-mismatches 2 --read-edit-dist 13 --max-insertion-length 20 --max-deletion-length 30 --num-threads 2")

command <- paste("tophat2 -o ", output.aligned, " --GTF ", genesets.loc, " --transcriptome-index ", static.dir," ", bowtie.index," ", sample.loc, " >& ", log.file)
if(!file.exists(output.aligned)) system(command, intern = TRUE)

### Sort the BAM (1 hour) ###
# -@ : Number of CPU cores to use
# -m : Memory per CPU core
outputSBAM     <- paste0(output.dir, sample.name, "_trim_align_sort.bam")

command        <- paste0(command, "samtools sort -@ 4 -m 2G -o ", output.aligned," ", sample.name, " > ", outputSBAM)
if(!file.exists(outputSBAM)) system(command, intern=TRUE)

### Index the sorted BAM file (10 minutes) ###
outputSBAI   <- paste0(output.dir, sample.name, "_trim_align_sort.bai")

command      <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
if(!file.exists(outputSBAI)) system(command, intern=TRUE)

### Mark duplicates, using the Picard tools
outputSBAID    <- paste0(output.dir, sample.name, "_trim_align_sort_dedup.bam")
outputMetrics  <- paste0(output.dir, sample.name, "_matrics.txt")

command <- paste0("java -jar MarkDuplicates.jar REMOVE_DUPLICATES=true INPUT=", outputSBAI, " OUTPUT=", outputSBAID," METRICS_FILE=", outputMetrics)
if(!file.exists(outputSBAID)) system(command, intern=TRUE)

### Get some basic statistics (5 to 10 minutes)
command      <- paste0("samtools flagstat ", outputSBAID)
system(command, intern=TRUE)


### TODO add the whole GATK part of the pipeline (see adjust the stuff below, so it matched the stuff above)
###- NEED TO CHECK FORM HERE ON

### Indel Realign
outputSNPS    <- "output.snp.intervals"
outputSIBAM   <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.bam")
knownindels   <- paste0(referenceDir, "/mgp.v3.indels.rsIDdbSNPv137.vcf")                              # Reference, download from: ftp://ftp-mouse.sanger.ac.uk/
if(!file.exists(outputSNPS)){
  command <- paste0("java -Xmx4g -jar ", gatk, " -nt 4 -T RealignerTargetCreator -R ", reference, " -known ", knownindels, " -o ", outputSNPS, " -U ALLOW_N_CIGAR_READS")
  cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n")   # Call the GATK RealignerTargetCreator, only need to be done because the knownSNPlist does not change
}
command <- paste0("java -Xmx4g -jar ", gatk, " -T IndelRealigner -R ", reference, " -targetIntervals ", outputSNPS, " -maxReads 150000 -I ", outputSBAID, " -o ",outputSIBAM," -known ",knownindels, " -U ALLOW_N_CIGAR_READS")
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n")     # Call the GATK IndelRealigner

### Base Recalibration
knownsnps     <- paste0(referenceDir, "/mgp.v3.snps.rsIDdbSNPv137.vcf")                                # Reference, download from: ftp://ftp-mouse.sanger.ac.uk/
covfile1      <- paste0(fileBase, "1.covariates")
covfile2      <- paste0(fileBase, "2.covariates")
plotfile      <- paste0(fileBase, "recalibration.pdf")
outputSIRBAM  <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam")

command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIBAM," -o ",  covfile1, " -U ALLOW_N_CIGAR_READS")
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n")   # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T PrintReads -R ", reference," -I ", outputSIBAM," -BQSR ", covfile1, " -U ALLOW_N_CIGAR_READS -o ", outputSIRBAM)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n")   # Call the GATK PrintReads
command <- paste0("java -Xmx4g -jar ", gatk, " -nct 4 -T BaseRecalibrator -R ", reference, " -knownSites ", knownsnps, " -I ", outputSIRBAM," -o ", covfile2, " -U ALLOW_N_CIGAR_READS")
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n")   # Call the GATK BaseRecalibrator
command <- paste0("java -Xmx4g -jar ", gatk, " -T AnalyzeCovariates -R ", reference, " -before ", covfile1, " -after ", covfile2, " -U ALLOW_N_CIGAR_READS -plots ", plotfile)
cat(system(command, intern=TRUE), file=logfile, append=TRUE, sep="\n")   # Call the GATK AnalyzeCovariates


