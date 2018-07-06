wd <- "/home/danny/NAS/Horse/RNA/RNAseq-Kabardian"
setwd(wd)

folders <- list.files(pattern="output")
files <- NULL
for(folder in folders){
  files <- c(files, paste0(wd, "/", folder, "/", list.files(folder, ".recal.bam")))
}

cat(paste0("nohup samtools mpileup -g -f /home/danny/genomes/EquCab2/Equus_caballus.EquCab2.dna.fa ", paste(files, collapse=" "), " -o population.bcf &"), "\n")
cat(paste0("nohup bcftools call -c population.bcf | ~/Github/bcftools/vcfutils.pl varFilter -d 10 - > population.vcf &"), "\n")

#nohup samtools mpileup -g -f /home/danny/genomes/EquCab2/Equus_caballus.EquCab2.dna.fa /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/11N.output/11N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/11V.output/11V.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/14N.output/14N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/14V.output/14V.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/15N.output/15N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/15V.output/15V.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/16N.output/16N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/16V.output/16V.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/21N.output/21N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/21V.output/21V.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/22N.output/22N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/22V.output/22V.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/23N.output/23N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/23V.output/23V.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/24N.output/24N.aln.sort.rmdup.rg.recal.bam /home/danny/NAS/Horse/RNA/RNAseq-Kabardian/24V.output/24V.aln.sort.rmdup.rg.recal.bam -o population.bcf &

#nohup /home/danny/Github/bcftools/bcftools call -c population.bcf | ~/Github/bcftools/vcfutils.pl varFilter -d 10 - > population.vcf &


### Gene expression
library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")

wd <- "/home/danny/NAS/Horse/RNA/RNAseq-Kabardian"
setwd(wd)

folders <- list.files(pattern="output")
files <- NULL
for(folder in folders){
  filename <- list.files(folder, ".recal.bam")
  if(length(filename) == 1) files <- c(files, paste0(wd, "/", folder, "/", filename))
}

# Create sizes.genome from fai index file:  cut -f1,2 ~/genomes/EquCab2/Equus_caballus.EquCab2.dna.fa.fai > sizes.genome
chrominfo <- read.table("/home/danny/RNAseqHorses/sizes.genome", sep="\t", header=FALSE, colClasses=c("character","integer"))
chrominfo <- cbind(chrominfo, FALSE)
colnames(chrominfo) <- c("chrom", "length", "is_circular")

horse <- makeTxDbFromGFF("/home/danny/genomes/EquCab2/Equus_caballus.EquCab2.90.gtf", format = "gtf", 
                                         organism="Equus caballus", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-90/gtf/equus_caballus/")
exonsByGene   <- exonsBy(horse, by = "gene")

bamfilesafter <- BamFileList(files, yieldSize = 1000000)
se            <- summarizeOverlaps(exonsByGene, bamfilesafter, mode="Union", singleEnd=TRUE, ignore.strand=TRUE, fragments=FALSE)
head(assay(se))                                                                             # Show the first 10 lines of the data matrix
write.table(assay(se), file="~/RNAseqHorses/rawreads.txt",sep="\t", quote=FALSE)

rawreads <- assay(se)
rawreads[rawreads == 0] <- NA                                                               # Change 0 reads to NA, so we can do quantile normalisation
rawreadsQnorm <- normalize.quantiles(as.matrix(rawreads))
rawreadsQnorm[is.na(rawreadsQnorm)] <- 0                                                    # Change back the NAs to 0 reads

# RPKM = (10^9 * C)/(N * L)
# C = Number of reads mapped to a gene
# N = Total mapped reads in the sample
# L = gene length in base-pairs for a gene

exonicGeneSizes <- lapply(exonsByGene, function(x){ sum(width(reduce(x))) })                # Get the length of each gene using only the exons
N <- apply(rawreadsQnorm, 2, sum)                                                           # Get the number of reads in all samples

n <- 1
RPKM <- t(apply(rawreadsQnorm, 1, function(C){                                              # Calculate the RPKM values per gene
  L     <- as.numeric(exonicGeneSizes[n])
  RPKM  <- (10^9 * C) / (N * L)
  n    <<- n + 1
  return(round(RPKM, d = 3))
}))

colnames(RPKM) <- colnames(rawreads)
rownames(RPKM) <- rownames(rawreads)
cat("We called expressions for", nrow(RPKM), "genes\n")

write.table(RPKM, file="~/RNAseqHorses/RPKMnormalized.txt", sep="\t")                                  # Write the normalized RPKM values to file
