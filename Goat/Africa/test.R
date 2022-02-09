setwd("D:/Edrive/Goat/")
sample.annot <- read.table(file = "Africa.Swiss.Adaptmap.samples.annotation.noDuplicates.txt", sep = "\t", quote="\"", colClasses = "character")

iix <- names(which(table(sample.annot[which(sample.annot[,"Continent"] == "Africa"),"Breed"]) > 35))


sample.annot[which(sample.annot[, "Breed"] %in% iix), "Country"]


which(sample.annot[, "Breed"] %in% iix & sample.annot[, "Country"] == "Burundi")