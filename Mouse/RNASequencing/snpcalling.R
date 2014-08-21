# snpcalling.R - Calling SNPs and indels using the GenomeAnalysisToolKit
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

# setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")

gatk          <- "/opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"
referenceDir  <- "genomes"
referenceName <- "Mus_musculus.GRCm38.74.dna"
reference     <- paste0(referenceDir, "/", referenceName, ".fasta")
knownsnps     <- paste0(referenceDir, "/", "mgp.v3.snps.rsIDdbSNPv137.vcf")                                # Reference, download from: ftp://ftp-mouse.sanger.ac.uk/

samples  <- read.table("FASTQ/sampledescription.txt", sep="\t", header=TRUE)
bamfiles <- dir()[grep("trimmed.aligned.sorted.realigned.recalibrated.bam", dir())]

for(bamfile in bamfiles){
  fileBase    <- strsplit(bamfiles[1],"_")[[1]][1]
  outputVCF   <- paste0(fileBase, ".snps.vcf")
  settings    <- "-recoverDanglingHeads -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0" 

  command <- paste0("java -Xmx4g -jar ", gatk, " -T HaplotypeCaller -R ", reference, " -I ", bamfile, "  --dbsnp ", knownsnps, " ", settings, " -o ", outputVCF)
  system(command, intern=TRUE)   # Call the GATK HaplotypeCaller
}
