#
# Wachow analysis Jan 2020
#

library(foreign)

setwd("D:/Edrive/Cow/Salma Wachow")
gts <- read.table("holstein_genotypes.txt", sep = "\t", header=TRUE, colClasses = "character")

setwd("D:/Edrive/Cow/HERDE/NeuDaten08/Herde/2019_11")
mydata <- read.dbf("MLP.DBF")
bestand <- read.dbf("BESTAND.DBF")

# Figure out MKG == NA or MKG == 0, and remove those from the data (no milk, no reliable SCS)
noData <- which(is.na(mydata[, "MKG"]) | mydata[, "MKG"] == 0)
mydata <- mydata[-noData,]

# Remove the NA cell counts
mydata <- mydata[-which(is.na(mydata[,"ZELLZAHL"])),]
mydata <- mydata[which(mydata[, "OHR"] %in% gts[, "Barcode.Ohmark"]),]

# How many genotyped cows do not have SCS data ?
length(which(!(gts[, "Barcode.Ohmark"] %in% mydata[, "OHR"])))
# Remove the cows without SCS
gts <- gts[-which(!(gts[, "Barcode.Ohmark"] %in% mydata[, "OHR"])),]

mmatrix <- c()
for(r in 1:nrow(mydata)){
  ohr <- mydata[r, "OHR"]
  inBestand <- which(as.character(bestand[, "OHR"]) == as.character(ohr))
  birthdate <- NA
  firstcalf <- NA
  if(length(inBestand) == 1){ 
    birthdate <- as.character(bestand[inBestand, "GEBURT"])
    firstcalf <- as.character(bestand[inBestand, "KALBUNG1"])
    firstcalfindays <- as.character(as.Date(firstcalf, format="%Y-%m-%d") - as.Date(birthdate, format="%Y-%m-%d"))
  }
  
  mlpdata <- c(as.character(mydata[r, "OHR"]), birthdate, firstcalfindays,
               as.character(mydata[r, "DATUM"]), 
               as.character(mydata[r, "LAKTATION"]), 
               as.character(mydata[r, "MKG"]),as.character(mydata[r, "FETT"]),as.character(mydata[r, "EIWEISS"]), as.character(mydata[r, "ZELLZAHL"]))
  mlpdata <- c(mlpdata, as.character(gts[which(as.character(gts[, "Barcode.Ohmark"]) == as.character(mydata[r, "OHR"])), -1]))

  mmatrix <- rbind(mmatrix, mlpdata)
}

colnames(mmatrix) <- c("ID", "Birthdate", "FirstCalfIndays","Sampledate", "Lactation", "MKG", "Fat", "Protein", "SCC", colnames(gts)[-1])
rownames(mmatrix) <- 1:nrow(mmatrix)

setwd("D:/Edrive/Cow/Salma Wachow")

write.table(mmatrix, "phenotypes2020.wachow.txt", sep = "\t", quote=FALSE)