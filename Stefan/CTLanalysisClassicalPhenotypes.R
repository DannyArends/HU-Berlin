# CTLanalysisClassicalPhenotypes.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Small script to perform CTL analysis on classical phenotypes

library(ctl)
setwd("d:/Stefan_Mouse_F3")

cross <- read.cross("csv", file="100419_c_f3_phegen.csv", sep=";")
geno <- pull.geno(cross)
pheno <- pull.pheno(cross)

traits <- c(4:6, 9:10, 12:17,22:26, 27,28,30:50, 53:54, 56:60)

pheno <- pheno[,traits]

derived <- grep("log", colnames(pheno))
derived <- c(derived, grep("drwgt_Rm", colnames(pheno)))
pheno <- pheno[,-derived]

ctlscan <- CTLscan(geno, pheno, verbose = TRUE)
significant <- ctl.lineplot(ctlscan)
