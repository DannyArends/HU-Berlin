# Analysis of Wellcome-CTC data from the mouse diversity array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#

setwd("E:/Mouse/DNA/DiversityArray/CTC")

snpdata1 <- read.table("RAW/Strains-0.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
snpdata2 <- read.table("RAW/Strains-200.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
snpdata3 <- read.table("RAW/Strains-400.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)

alldata <- cbind(snpdata1, snpdata2, snpdata3)[, -c(203, 204, 405, 406)]
write.table(alldata, file="Analysis/alldata.txt", sep="\t")