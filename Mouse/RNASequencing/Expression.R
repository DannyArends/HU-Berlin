# RNA Seq - Expression data analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

# Download GTF file from: ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/

library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")

setwd("E:/Mouse/RNA/Sequencing")
mouse         <- makeTranscriptDbFromGFF("GTF/Mus_musculus.GRCm38.76.gtf", format = "gtf")
exonsByGene   <- exonsBy(mouse, by = "gene")
infiles       <- list.files("Analysis", pattern="recalibrated.bam$", full=TRUE)
bamfiles      <- BamFileList(infiles, yieldSize=100000)
se            <- summarizeOverlaps(exonsByGene, bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

head(assay(se))

#C = Number of reads mapped to a gene
#N = Total mapped reads in the experiment
#L = gene length in base-pairs for a gene

exonicGeneSizes <- lapply(exonsByGene, function(x){sum(width(reduce(x)))})
N <- apply(assay(se), 2, sum)
RPKM = (10^9 * C)/(N * L)
