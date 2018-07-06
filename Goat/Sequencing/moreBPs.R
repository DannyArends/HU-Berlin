setwd("D:/Edrive/Goat/DNA/Sequencing")
library(seqinr)
library(Biostrings)

fastaSeq <- read.fasta("genome/GoatGenome_6_8_14.fa")

bparound <- NULL

sequence <- fastaSeq[["ENA|CM004567|CM004567.1"]][(85981710 - 400) : (85981710 + 400)]
around <- paste0(paste0(sequence[1:400],collapse=""), "[", sequence[401], "]", paste0(sequence[402:801],collapse=""))
bparound <- rbind(bparound, c("ENA|CM004567|CM004567.1",85981710, sequence[401], around))

sequence <- fastaSeq[["ENA|CM004567|CM004567.1"]][(85982615 - 400) : (85982615 + 400)]
around <- paste0(paste0(sequence[1:400],collapse=""), "[", sequence[401], "]", paste0(sequence[402:801],collapse=""))
bparound <- rbind(bparound, c("ENA|CM004567|CM004567.1", 85982615, sequence[401], around))

sequence <- fastaSeq[["ENA|CM004567|CM004567.1"]][(86008016 - 400) : (86008016 + 400)]
around <- paste0(paste0(sequence[1:400],collapse=""), "[", sequence[401], "]", paste0(sequence[402:801],collapse=""))
bparound <- rbind(bparound, c("ENA|CM004567|CM004567.1", 86008016, sequence[401], around))

sequence <- fastaSeq[["ENA|CM004567|CM004567.1"]][(86079098 - 400) : (86079098 + 400)]
around <- paste0(paste0(sequence[1:400],collapse=""), "[", sequence[401], "]", paste0(sequence[402:801],collapse=""))
bparound <- rbind(bparound, c("ENA|CM004567|CM004567.1", 86079098, sequence[401], around))

sequence <- fastaSeq[["ENA|CM004567|CM004567.1"]][(86208960 - 400) : (86208960 + 400)]
around <- paste0(paste0(sequence[1:400],collapse=""), "[", sequence[401], "]", paste0(sequence[402:801],collapse=""))
bparound <- rbind(bparound, c("ENA|CM004567|CM004567.1", 86208960, sequence[401], around))

sequence <- fastaSeq[["ENA|CM004575|CM004575.1"]][(81331681 - 400) : (81331681 + 400)]
around <- paste0(paste0(sequence[1:400],collapse=""), "[", sequence[401], "]", paste0(sequence[402:801],collapse=""))
bparound <- rbind(bparound, c("ENA|CM004567|CM004575.1", 81331681, sequence[401], around))

sequence <- fastaSeq[["ENA|CM004575|CM004575.1"]][(81326608 - 400) : (81326608 + 400)]
around <- paste0(paste0(sequence[1:400],collapse=""), "[", sequence[401], "]", paste0(sequence[402:801],collapse=""))
bparound <- rbind(bparound, c("ENA|CM004567|CM004575.1", 81326608, sequence[401], around))


write.table(bparound, "400bparound.txt", sep="\t", row.names=F, col.names=FALSE, quote=F)