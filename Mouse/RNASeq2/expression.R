
# Loading data
setwd("E:/Mouse/RNA/ArrayDesign/Atlas data")
arrays <- read.table("Annotation/arrays.txt", header=TRUE, sep="\t", colClasses="character")
alldata <- read.table("Analysis/geneexpression.txt", sep="\t", header=TRUE)