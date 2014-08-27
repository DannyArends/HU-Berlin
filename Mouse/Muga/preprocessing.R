# Preprocessing of the MegaMuga data, recreate the genotype call matrix from the raw read data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/DNA/MegaMuga/Humboldt_Berlin_MEGMUGV01_20140817")

con  <- file("Humboldt_Berlin_MEGMUGV01_20140817_FinalReport.txt", open = "r")

fileheader   <- readLines(con, n = 9, warn = FALSE)
matrixheader <-  unlist(strsplit(readLines(con, n = 1, warn = FALSE), "\t"))

DATA <- NULL
cnt <- 1
while(length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0 && cnt < 100){
  DATA <- rbind(DATA, unlist(strsplit(oneLine,"\t")))
  cnt <- cnt + 1
}
close(con)
colnames(DATA) <- matrixheader

