# RNA Seq - Expression data analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015

library(Rsamtools)

bamcoverage <- function (bamfile) {
  param <- ScanBamParam(what=c("pos","qwidth","qname","strand"))

  bam <- scanBam(bamfile, maxMemory=1024*3, param=param)[[1]]                                                    # Read in the bam file, the result comes in nested lists
  ind <- !is.na(bam$pos)                                                          # Filter reads without match position

  bam <- lapply(bam, function(x) x[ind])                                          # Remove non-matches, they are not relevant to us
  ranges <- IRanges(start=bam$pos, width=bam$qwidth, names=make.names(bam$qname, unique=TRUE))
  ranges <- GRanges(seqnames=Rle(bam$rname), ranges=ranges, strand=Rle(bam$strand), flag=bam$flag, readid=bam$rname)

  return (mean(coverage(ranges)))      
}

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI")

bamfiles <- list.files("Analysis","*recalibrated.bam$")
bamcoverage(paste0("Analysis/", bamfiles[1]))