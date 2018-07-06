# Pipeline for single end RNA-Seq analysis
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
readgroupID  <- as.numeric(gsub("V","",gsub("N", "", fname)))

cat(base,", ", fname, ", ", filebase, ", ", readgroupID,"\n")

# Locations, executable names and constants
reference              <- "/home/danny/genomes/EquCab2/Equus_caballus.EquCab2.dna"
reference.fa           <- paste0(reference, ".fa")        #!
reference.fa.gz        <- paste0(reference, ".fa.gz")     #!
reference.fa.fai       <- paste0(reference, ".fa.fai")
reference.dict         <- paste0(reference, ".dict")
reference.bt2idx       <- paste0(reference, "")
reference.bt2idx.fa    <- paste0(reference, ".fa") #!

reference.snps      <- "/home/danny/genomes/EquCab2/Equus_caballus.vcf.gz"
reference.gtf       <- "/home/danny/genomes/EquCab2/Equus_caballus.EquCab2.90.gtf.gz"

picarddict.exec     <- "java -Xmx4g -jar /home/danny/picard-tools-1.119/CreateSequenceDictionary.jar"
picardrmdup.exec    <- "java -Xmx4g -jar /home/danny/picard-tools-1.119/MarkDuplicates.jar"
picardrg.exec       <- "java -Xmx4g -jar /home/danny/picard-tools-1.119/AddOrReplaceReadGroups.jar"
picardreo.exec      <- "java -Xmx4g -jar /home/danny/picard-tools-1.119/ReorderSam.jar"
gatk.exec           <- "java -Xmx16g -jar /home/danny/Github/GATK-3.7/GenomeAnalysisTK.jar"
b2build.exec        <- "/home/danny/bowtie2-2.3.2/bowtie2-build"
bowtie2.exec        <- "/home/danny/bowtie2-2.3.2/bowtie2"
samtools.exec       <- "samtools"
tabix.exec          <- "tabix"

transcriptome.index  <- "/home/danny/genomes/EquCab2/tophat2/EquCab2.genes"

# Trim reads by Trimmomatic
trimmomatic.files  <- c(paste0(inputfile,"_R1.fastq.gz"), paste0(filebase,"_R1.P.fastq.gz"))
trimmomatic.exec   <- "java -jar /home/danny/Github/trimmomatic/trimmomatic-0.36/dist/jar/trimmomatic-0.36.jar"
trimmomatic.opts   <- "ILLUMINACLIP:/home/danny/Github/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

execute(paste0(trimmomatic.exec, " SE ", paste0(trimmomatic.files, collapse=" "), " ", trimmomatic.opts), trimmomatic.files[2])

# Aligment using tophat2
alignment.file  <- paste0(filebase, ".tophat2.aln.bam")
tophat.options  <- paste0("--rg-id ", readgroupID, " --rg-sample ", fname, " --rg-library RNA-seq --rg-platform Illumina --b2-sensitive --num-threads 4", collapse="")
tophat.log      <- paste0(filebase, ".tophat2.log")
execute(paste0("/home/danny/Github/tophat-2.1.1.Linux_x86_64/tophat2"," -o ", alignment.file, " --transcriptome-index ", transcriptome.index, " ", tophat.options, " ", reference.bt2idx, " ", trimmomatic.files[2], " > ", tophat.log, collapse=""), alignment.file)
alignment.file    <- paste0(filebase, ".tophat2.aln.bam/accepted_hits.bam")

# Sort and remove duplicates
samtools.file    <- paste0(filebase, ".aln.sort")
samtools.bamfile <- paste0(samtools.file, ".bam")
picard.file <- paste0(filebase, ".aln.sort.rmdup.bam")
metrics.file <- paste0(filebase, ".metrics.txt")

execute(paste0(samtools.exec, " sort -@ 4 -m 2G ",alignment.file, " -o ", samtools.bamfile), samtools.bamfile)
execute(paste0(picardrmdup.exec, " REMOVE_DUPLICATES=true INPUT=", samtools.bamfile, " OUTPUT=", picard.file, " METRICS_FILE=", metrics.file), picard.file)
#execute(paste0(samtools.exec, " flagstat ", picard.file))

# Add read groups
picard.rgfile             <-  paste0(filebase, ".aln.sort.rmdup.rg.bam")
picard.rgfile.reordered   <-  paste0(filebase, ".aln.sort.rmdup.rg.o.bam")
execute(paste0(picardrg.exec, " INPUT=", picard.file, " OUTPUT=", picard.rgfile, " CREATE_INDEX=false RGID=",readgroupID," RGLB=LIB",readgroupID," RGPL=Illumina RGPU=X RGSM=",readgroupID), picard.rgfile)

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
execute(paste0(gatk.exec, " -nct 4 -T BaseRecalibrator -R ", reference.fa," -knownSites ", reference.snps, " -I ",picard.rgfile.reordered, " -o ", gatk.cov1," -U ALLOW_N_CIGAR_READS"), gatk.cov1)
gatk.recal <- paste0(filebase, ".aln.sort.rmdup.rg.recal.bam")
execute(paste0(gatk.exec, " -nct 4 -T PrintReads -R ", reference.fa," -I ", picard.rgfile.reordered," -BQSR ", gatk.cov1," -U ALLOW_N_CIGAR_READS -o ", gatk.recal), gatk.recal)
gatk.cov2 <- paste0(filebase, ".2.covariates")
execute(paste0(gatk.exec, " -nct 4 -T BaseRecalibrator -R ", reference.fa," -knownSites ", reference.snps, "  -I ",gatk.recal," -o ", gatk.cov2," -U ALLOW_N_CIGAR_READS"), gatk.cov2)
gatk.plots <- paste0(filebase, ".plots.pdf")
execute(paste0(gatk.exec, " -T AnalyzeCovariates -R ", reference.fa," -before ", gatk.cov1, " -after ", gatk.cov2," -U ALLOW_N_CIGAR_READS -plots ", gatk.plots), gatk.plots)

# Done with the whole RNA-Seq pipeline
q("no")
