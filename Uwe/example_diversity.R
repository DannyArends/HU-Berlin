# Example of using diversity to calculate some basic and more advanced popualtion statistics
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends

setwd("D:/Github/HU-Berlin/Uwe")

install.packages("diveRsity")               # This is only required once, after installing you do not have to install again
library(diveRsity)

basicStats <- divBasic(infile = "genotypes_genpop.txt", outfile="fstOnlyOut.txt", gp=2, bootstraps = 100)
names(basicStats$fis) <- colnames(basicStats$Ho)

basicStats$Ho["overall",]                   # Observed
basicStats$He["overall",]                   # Expected
basicStats$fis[["Tagg,"]]["overall",]       # Fis
basicStats$fis[["Ni,"]]["overall",]         # Fis
basicStats$fis[["Nu,"]]["overall",]         # Fis
basicStats$fis[["Dese,"]]["overall",]       # Fis
basicStats$fis[["Combined"]]["overall",]    # Fis

advancedStats <- diffCalc(infile = "genotypes_genpop.txt", outfile="fstOnlyOut.txt", fst=TRUE, pairwise=TRUE, boots = 1000)

advancedStats$pairwise$Fst                  # Pairwise Fst per populations

