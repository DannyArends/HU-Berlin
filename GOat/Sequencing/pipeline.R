# Pipeline for DNA re-seq analysis on chicken
#
# copyright (c) 2015-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

#cmdlineargs   <- commandArgs(trailingOnly = TRUE)

samples <- read.table("sample_names.txt",header=TRUE)
reference   <- "genome/GoatGenome_6_8_14.fa"

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

for(x in 1:nrow(samples)) {
  fileBase    <- as.character(samples[x, 1])
  inputfiles  <- c(paste0("fastq/DeMux/", fileBase, "_R1.fastq"), paste0("fastq/DeMux/", fileBase, "_R2.fastq"))

  ########### ANALYSIS ###########

  ### BWA: Alignment against genome pipe the output to the next step (SAM to BAM)
  # -v 2  : Verbose level 2 only errors and warnings
  # -t 6  : Number of threads to use with BWA
  command     <- paste0("/home/danny/Github/bwa/bwa-0.7.15/bwa mem -v 2 -t 6 -A 3 -B 2 -U 4 -O 2 -E 0 -T 10 -a ", reference," ", inputfiles[1], " ", inputfiles[2], " | samtools view -Sb - > alignments/", fileBase, ".bam")
  execute(command)

  ### Sort the aligned BAM file
  outputSBAM <- paste0("alignments/", fileBase, ".sorted.bam")
  command    <- paste0("samtools sort -@ 4 -m 2G alignments/", fileBase, ".bam > ", outputSBAM)
  execute(command)

  ### Index the aligned BAM file
  outputSBAI <- paste0("alignments/", fileBase, ".sorted.bai")
  command    <- paste0("samtools index ", outputSBAM, " ", outputSBAI)
  execute(command)

  ### Get some basic statistics regarding the alignment
  command      <- paste0("samtools flagstat ", outputSBAM)
  execute(command)
}

q("no")

### Call SNPs use reads with a minimum mapping quality of 10, and use the multiallelic caller from BCFtools, then filter using varFilter

#ENA|CM004567|CM004567.1	85973463	85996270	+	CSN1S1
#ENA|CM004567|CM004567.1:85973463-85996270
samtools mpileup -q 10 -r "ENA|CM004567|CM004567.1:85973463-85996270" -t DP -uf genome/GoatGenome_6_8_14.fa alignments/*.rg.bam | /home/danny/Github/bcftools/bcftools call -m - > CSN1S1.raw.bcf
/home/danny/Github/bcftools/bcftools view CSN1S1.raw.bcf | /home/danny/Github/bcftools/vcfutils.pl varFilter  > CSN1S1.SNPs.vcf

#ENA|CM004567|CM004567.1	86001250	86016321	+	CSN2
#ENA|CM004567|CM004567.1:86001250-86016321
samtools mpileup -q 10 -r "ENA|CM004567|CM004567.1:86001250-86016321" -t DP -uf genome/GoatGenome_6_8_14.fa alignments/*.rg.bam | /home/danny/Github/bcftools/bcftools call -m - > CSN2.raw.bcf
/home/danny/Github/bcftools/bcftools view CSN2.raw.bcf | /home/danny/Github/bcftools/vcfutils.pl varFilter  > CSN2.SNPs.vcf

#ENA|CM004567|CM004567.1	86071845	86094539	+	CSN1S2
#ENA|CM004567|CM004567.1:86071845-86094539
samtools mpileup -q 10 -r "ENA|CM004567|CM004567.1:86071845-86094539" -t DP -uf genome/GoatGenome_6_8_14.fa alignments/*.rg.bam | /home/danny/Github/bcftools/bcftools call -m - > CSN1S2.raw.bcf
/home/danny/Github/bcftools/bcftools view CSN1S2.raw.bcf | /home/danny/Github/bcftools/vcfutils.pl varFilter  > CSN1S2.SNPs.vcf

#ENA|CM004567|CM004567.1	86192263	86212376	+	CSN3
#ENA|CM004567|CM004567.1:86192263-86212376
samtools mpileup -q 10 -r "ENA|CM004567|CM004567.1:86192263-86212376" -t DP -uf genome/GoatGenome_6_8_14.fa alignments/*.rg.bam | /home/danny/Github/bcftools/bcftools call -m - > CSN3.raw.bcf
/home/danny/Github/bcftools/bcftools view CSN3.raw.bcf | /home/danny/Github/bcftools/vcfutils.pl varFilter  > CSN3.SNPs.vcf

#ENA|CM004569|CM004569.1	70790018	70807076	+	STC1
#ENA|CM004569|CM004569.1:70790018-70807076
samtools mpileup -q 10 -r "ENA|CM004569|CM004569.1:70790018-70807076" -t DP -uf genome/GoatGenome_6_8_14.fa alignments/*.rg.bam | /home/danny/Github/bcftools/bcftools call -m - > STC1.raw.bcf
/home/danny/Github/bcftools/bcftools view STC1.raw.bcf | /home/danny/Github/bcftools/vcfutils.pl varFilter  > STC1.SNPs.vcf

#ENA|CM004575|CM004575.1	81324989	81339811	+	DGAT1
#ENA|CM004575|CM004575.1:81324989-81339811
samtools mpileup -q 10 -r "ENA|CM004575|CM004575.1:81324989-81339811" -t DP -uf genome/GoatGenome_6_8_14.fa alignments/*.rg.bam | /home/danny/Github/bcftools/bcftools call -m - > DGAT1.raw.bcf
/home/danny/Github/bcftools/bcftools view DGAT1.raw.bcf | /home/danny/Github/bcftools/vcfutils.pl varFilter  > DGAT1.SNPs.vcf
