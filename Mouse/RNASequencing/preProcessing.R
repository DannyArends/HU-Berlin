# Preprocessing for RNA seq analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

########### SETUP THE REFERENCE AND KNOWN SNPS AND INDELS ###########

referenceDir  <- "genomes"
referenceName <- "Mus_musculus.GRCm38.74.dna"

### Download the reference genome, and create a single big FASTA
reference     <- paste0(referenceDir, "/", referenceName, ".fasta")
chromosomes <- c(1:19, "X", "Y", "MT")
cat("", file=reference)
for(chr in chromosomes){
  system(paste0("wget -P temp ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa.gz"))
  system(paste0("gunzip temp/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa.gz"))
  fastadata <- readLines(paste0("temp/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa"))
  cat(fastadata, sep="\n", file = reference, append = TRUE)
  system(paste0("rm temp/Mus_musculus.GRCm38.74.dna.chromosome.", chr, ".fa.gz")) # Delete the temp fasta file
}

### Download the known INDELS and SNPs
system(paste0("wget -P ", referenceDir, " ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/mgp.v3.indels.rsIDdbSNPv137.vcf.gz"))
system(paste0("wget -P ", referenceDir, " ftp://ftp-mouse.sanger.ac.uk/REL-1303-SNPs_Indels-GRCm38/mgp.v3.snps.rsIDdbSNPv137.vcf.gz"))

########### PRE-PROCESSING THE REFERENCE ###########

### Create A BWA index file from the reference genome / transcriptome (10 to 100 minutes) ###
system(paste0("/opt/bwa-0.7.10/bwa index ", reference))

### Create A SAMTOOLS index file from the reference genome ###
system(paste0("samtools faidx ", reference))

### Create a Picard dictionary file
referenceDict <- paste0(referenceDir, "/", referenceName, ".dict")
system(paste0("java -jar /opt/picard-tools-1.99/CreateSequenceDictionary.jar R=", reference, " O=", referenceDict))
