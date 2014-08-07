# Use the ShortRead library to create subsets of the big Fastq files
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

library(ShortRead)
  
setwd("E:/Mouse/RNA/Sequencing/FASTQ")
f1 <- FastqSampler("4868_GCCAAT_L001_R1_001.fastq.gz", n=1000000)
f2 <- FastqSampler("4868_GCCAAT_L001_R2_001.fastq.gz", n=1000000)

set.seed(123L); p1 <- yield(f1)
set.seed(123L); p2 <- yield(f2)

writeFastq(p1, "4868_GCCAAT_L001_R1_001_S.fastq.gz")
writeFastq(p2, "4868_GCCAAT_L001_R2_001_S.fastq.gz")