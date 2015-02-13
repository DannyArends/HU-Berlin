library("GenomicAlignments")
library("GenomicFeatures")
library("Rsamtools")
library("preprocessCore")
library("BiocParallel")

# Chromosome info structure, holding the number of chromosomes, their names and the length per chromosome
chrominfo     <- read.table("Mouse/GTF/MouseChrInfo.txt", sep="\t", header=TRUE, colClasses=c("character","integer","logical"))

# Download GTF file from: ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/
mouse         <- makeTranscriptDbFromGFF("Mouse/GTF/Mus_musculus.GRCm38.76.gtf", format = "gtf", exonRankAttributeName="exon_number", 
                                         species="Mus musculus", chrominfo=chrominfo, dataSource="ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/")

transcriptsByGene <- transcriptsBy(mouse, by = "exon")
infiles           <- list.files(path = "./Mouse", pattern="recalibrated.bam$", full=TRUE)
bamfiles          <- BamFileList(infiles, yieldSize = 1000000, asMates=TRUE)

register(MulticoreParam(workers = 5))

se                <- summarizeOverlaps(transcriptsByGene, bamfiles, mode="Union", singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE)

head(assay(se))                                                                             # Show the first 10 lines of the data matrix
rawreads <- assay(se)                                                                       # Extract the raw-reads per gene

write.table(rawreads, file=paste0("Mouse/RawReadsPerExon.txt"), sep="\t")                         # Save the raw reads to a file

q("no")