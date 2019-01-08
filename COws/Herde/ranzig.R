# Herde matching script Junling
#
library(foreign)
setwd("D:/Edrive/Cow/HERDE/NeuDaten08/")

ranzigLKV <- read.table("TabellePatricia/Ranzig.txt", sep = "\t", header=TRUE)
ranzigLKV[,"Datum"] <- as.Date(ranzigLKV[,"Datum"], "%m/%d/%Y")

D2018_2 <- cbind(F = "2018_2", read.dbf("Herde/2018_02/BESTAND.DBF"))
D2018_7 <- cbind(F = "2018_7", read.dbf("Herde/2018_07/BESTAND.DBF"))
D2018_L  <- cbind(F = "2018_L", read.dbf("Herde/HW_91237082_180830/BESTAND.DBF"))

#allD <- rbind(D2018_2, D2018_7)
allD <- rbind(D2018_L)
allD <- allD[, c("F", "OHR", "STALL", "GEBURT", "KALBUNG1", "ABGANG")]

duplicatedEntries <- which(duplicated(allD))
if(length(duplicatedEntries) > 0) allD <- allD[-duplicatedEntries,]

allD[,"GEBURT"] <- as.Date(allD[,"GEBURT"], "%Y/%m/%d")
allD[,"KALBUNG1"] <- as.Date(allD[,"KALBUNG1"], "%Y/%m/%d")
allD[,"ABGANG"] <- as.Date(allD[,"ABGANG"], "%Y/%m/%d")

good <- 0
bad <- 0
miss <- 0
for(x in 1:nrow(ranzigLKV)){
  LKV_Date <- ranzigLKV[x, "Datum"]
  sameStall <- which(allD[,"STALL"] == ranzigLKV[x, "Stallnummer"])
  allD[sameStall, ]
  firstkalb <- allD[sameStall, "KALBUNG1"]
  abgang <- allD[sameStall, "ABGANG"]
  abgang[is.na(abgang)] <- Sys.Date()  # Set the Abgang to today if the cow is still here
  ii <- which(firstkalb < LKV_Date & LKV_Date < abgang)
  estCowID <- as.character(unique(allD[sameStall, ][ii, "OHR"]))
  cat(x, " ", ranzigLKV[x, "Stallnummer"], " ", as.character(LKV_Date), " ", estCowID , "\n")
  if(length(estCowID) == 1){
    ranzigLKV[x,"Ohrnummer"] <- estCowID
    good <- good + 1
  }
  if(length(estCowID) > 1) bad <- bad + 1
  if(length(estCowID) == 0) miss <- miss + 1
}

piedata <- c(good / nrow(ranzigLKV),bad / nrow(ranzigLKV), miss / nrow(ranzigLKV))
perc <- round(piedata * 100,1)
names(piedata) <- paste0(c("Matched ", "Multiple matches ", "No match "), perc, "%")
pie(piedata, col=c("green", "gray", "darkgray"), main = "Cows identification between LKV and HERDE\nBetrieb 1")

K2018_2 <- read.dbf("Herde/2018_02/KRANK.DBF")
K2018_7b <- read.dbf("Herde/2018_07b/KRANK.DBF")

allK <- rbind(K2018_2, K2018_7b)
allK <- allK[, c("OHR", "DATUM", "STAUFEN")]

duplicatedEntries <- which(duplicated(allK))
if(length(duplicatedEntries) > 0) allK <- allK[-duplicatedEntries,]
allK <- allK[grep("1.13", allK[,"STAUFEN"]),]

ranzigLKV <- cbind(ranzigLKV, staufen = NA)

LKVandHerde <- 0
onlyLKV <- 0
notidentified <- 0
matches <- c()
for (x in 1:nrow(ranzigLKV)) {
  LKV_Date <- ranzigLKV[x, "Datum"]
  ohr <- ranzigLKV[x, "Ohrnummer"]
  if(is.na(ohr)){
    notidentified <- notidentified + 1
  }else{
    sameOhr <- allK[which(allK[,"OHR"] == ohr),]
    from <- LKV_Date - 14
    to <- LKV_Date + 14
    ii <- which(from < sameOhr[,"DATUM"] & sameOhr[,"DATUM"] < to)
    if(length(ii) > 0){
      staufen <- as.character(unique(sameOhr[ii, "STAUFEN"]))
      cat(x, " mastites in HERDE,", staufen,"\n")
      ranzigLKV[x, "staufen"] <- paste0(staufen, collapse=", ")
      LKVandHerde <- LKVandHerde + 1
      matches <- c(matches, x)
    }else{
      cat(x, " no mastites in HERDE\n")
      onlyLKV <- onlyLKV + 1
    }
  }
}

piedata <- c(LKVandHerde / nrow(ranzigLKV), onlyLKV / nrow(ranzigLKV), notidentified / nrow(ranzigLKV))
perc <- round(piedata * 100,1)
names(piedata) <- paste0(c("Has Entry ", "No entry in HERDE ", "No / Multiple match "), perc, "%")
pie(piedata, col=c("green", "red", "gray"), main = "Mastites between LKV and HERDE\nBetrieb 1")

multipleStaufen <- length(which(lapply(strsplit(ranzigLKV[matches, "staufen"],", "), length) > 1))
oneStaufen <- which(lapply(strsplit(ranzigLKV[matches, "staufen"],", "), length) == 1)

staufenTable <- table(ranzigLKV[matches, "staufen"][oneStaufen])

piedata <- c(round(staufenTable / nrow(ranzigLKV) * 100 ,1), 
  "Multiple staufen" = round(multipleStaufen / nrow(ranzigLKV) * 100 ,1),
  "No entry in HERDE " = round((onlyLKV) / nrow(ranzigLKV) * 100 ,1),
  "No / Multiple match " = round((notidentified ) / nrow(ranzigLKV) * 100 ,1)
  )
perc <- round(piedata, 1)
names(piedata) <- paste0(names(piedata)," (", perc, "%)")
colz <- c(colorRampPalette(c("lightgreen", "darkgreen"))(length(staufenTable)), "purple", "red", "gray")
pie(piedata, col = colz, main = "Staufen in Herde\nBetrieb 1")

