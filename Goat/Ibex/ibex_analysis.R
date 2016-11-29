# Analysis of the Ibes Goat SNP data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2016
# first written Nov, 2016

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
chir1 <- read.csv("FilteredLocationCHIR1.0.txt", sep="\t", row.names=1)
chir1 <- cbind(chir1, Pos = (chir1[,"Start"] + chir1[,"Stop"])/2)

map <- chir1[,c("chrN", "Pos")]

setwd("D:/Edrive/Goat/DNA/Ibex")
snpdata <- read.csv("Goat_Nov2016_FinalReport.txt", sep='\t', skip=9, header=TRUE, na.strings=c("","--", "na", "NA"), row.names=1, colClasses="character")
colnames(snpdata) <- paste0("Ind", 1:48)

samples <- read.csv("Goat_Nov2016_Samples.txt", sep="\t", row.names = 1, colClasses="character")
rownames(samples) <- paste0("Ind", 1:48)

snpinfo <- read.csv("Goat_Nov2016_SNP.txt",sep="\t", row.names=1, colClasses="character")

notOnMap <- rownames(snpdata)[which(!(rownames(snpdata) %in% rownames(map)))]

snpdata <- snpdata[-which(rownames(snpdata) %in% notOnMap), ]
snpinfo <- snpinfo[-which(rownames(snpinfo) %in% notOnMap), ]

notOnArray <- rownames(map)[which(!rownames(map) %in% rownames(snpdata))]

map <- map[with(map, order(chrN, Pos)), ]
snpdata <- snpdata[rownames(map), ]
snpinfo <- snpinfo[rownames(map), ]

dim(map); map[1:5,]
dim(snpdata); snpdata[1:5,1:5]
dim(snpinfo); snpinfo[1:5,]

tooMuchMissing <- which(apply(apply(snpdata,1, is.na),2,sum) / ncol(snpdata) > 0.75)
tooLowGenTrain <- rownames(snpinfo)[which(snpinfo[,"GenTrain.Score"] < 0.6)]
tooLowMinorA <- rownames(snpinfo)[which(snpinfo[,"Minor.Freq"] < 0.05)]

badMarkers <- c(tooMuchMissing, tooLowGenTrain, tooLowMinorA)

map <-  map[-which(rownames(map) %in% badMarkers), ]
snpdata <- snpdata[-which(rownames(snpdata) %in% badMarkers), ]
snpinfo <- snpinfo[-which(rownames(snpinfo) %in% badMarkers), ]

dim(map); map[1:5,]
dim(snpdata); snpdata[1:5,1:5]
dim(snpinfo); snpinfo[1:5,]

snpinfo <- cbind(snpinfo, minorallele = NA, majorallele = NA)
for(x in rownames(snpdata)){
  rr <- table(unlist(strsplit(unlist(snpdata[x,]), "")))
  snpinfo[x, "minorallele"] <- names(which.min(rr))
  snpinfo[x, "majorallele"] <- names(which.max(rr))
}

numsnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames=list(rownames(snpdata), colnames(snpdata)))
for(x in rownames(snpdata)){
  minor.Homo <- paste0(snpinfo[x, "minorallele"], snpinfo[x, "minorallele"])
  major.Homo <- paste0(snpinfo[x, "majorallele"], snpinfo[x, "majorallele"])
  numsnpdata[x, snpdata[x,] == minor.Homo] <- 1
  numsnpdata[x, snpdata[x,] == major.Homo] <- 3
  numsnpdata[x, snpdata[x,] != major.Homo & snpdata[x,] != minor.Homo & !is.na(snpdata[x,])] <- 2
}

colnames(numsnpdata) <- samples[colnames(snpdata),"Origin.Species"]
distances <- dist(t(numsnpdata))
plot(hclust(distances))

absnpdata <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames=list(rownames(snpdata), colnames(snpdata)))
for(x in rownames(snpdata)){
  minor.Homo <- paste0(snpinfo[x, "minorallele"], snpinfo[x, "minorallele"])
  major.Homo <- paste0(snpinfo[x, "majorallele"], snpinfo[x, "majorallele"])
  absnpdata[x, snpdata[x,] == minor.Homo] <- "AA"
  absnpdata[x, snpdata[x,] == major.Homo] <- "BB"
  absnpdata[x, snpdata[x,] != major.Homo & snpdata[x,] != minor.Homo & !is.na(snpdata[x,])] <- "AB"
}

stammpinput <- t(absnpdata)
stammpinput <- cbind(Sample = rownames(stammpinput), Pop = as.character(samples[rownames(stammpinput), "Short"]), Ploidy = 2, Format = "BiA", stammpinput)
stammpinput <- as.data.frame(stammpinput)

stammpinput.freq <- stamppConvert(stammpinput, "r") # Frequencies
stammp.D.pop <- stamppNeisD(stammpinput.freq, TRUE) # Population D values
stammp.D.ind <- stamppNeisD(stammpinput.freq, FALSE) # Population D values
stammpinput.fst <- stamppFst(stammpinput.freq, 1000, 95, 4) # Population Fst values
stammpinput.fst$Fsts
write.table(stammpinput.fst$Fsts, file = "fsts.txt", sep = "\t")
