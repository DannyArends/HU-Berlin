# Herde matching script Junling
#

library(foreign)
setwd("D:/Edrive/Cow/HERDE/NeuDaten08/")

#Bestand.dbf file contains all animals day of arrival (GEBURT) and date of slaughter (ABGANG)
ranzig.bestand <- read.dbf("Herde/HW_91237082_180830/BESTAND.DBF")
ranzig.bestand <- ranzig.bestand[, c("OHR", "STALL", "GEBURT", "ABGANG")]

# Change to R-date object
ranzig.bestand[,"GEBURT"] <- as.Date(ranzig.bestand[,"GEBURT"], "%Y/%m/%d")
ranzig.bestand[,"ABGANG"] <- as.Date(ranzig.bestand[,"ABGANG"], "%Y/%m/%d")

# Diseases in the ranzig herde
ranzig.krank <- read.dbf("Herde/HW_91237082_180830/krank.DBF")

# Need a list of samples in the freezer at the HU
probes <- read.table("Herde/probesAtHU.txt", sep='\t', header=TRUE, colClasses="character")
probes[, "ID"] <- gsub("DE", "DE00", probes[, "ID"])
ranzig.probes <- probes[probes[, "Farm"] == "Ranzig",]
ranzig.probes <- ranzig.probes[which(ranzig.probes[,"ID"] != ""),]

ranzig.bestand <- ranzig.bestand[which(ranzig.bestand[, "OHR"] %in% probes[,"ID"]), ]
ranzig.probes <- ranzig.probes[which(ranzig.probes[,"ID"] %in% ranzig.bestand[, "OHR"]), ]

# For some reason there are duplicated samples of the same cow, we don't want to deal with them so we just remove them
ranzig.probes <- ranzig.probes[-which(duplicated(ranzig.probes[,"ID"])),]
# Take only the columns we are interested in 
ranzig.probes <- ranzig.probes[, c("ID", "TubeNr", "Code", "Farm", "Arrival")]
dim(ranzig.probes)

# Match the krank data to the tubes at the HU
ranzig.krank <- ranzig.krank[which(ranzig.krank[, "OHR"] %in% probes[,"ID"]), ]

# Figure out the mastites from Herde
ranzig.krank <- ranzig.krank[grep("1.13", ranzig.krank[,"STAUFEN"]),]
dim(ranzig.krank)

ranzig.probes <- cbind(ranzig.probes, MastitisInHerde = NA, Lactation = NA)
for(x in 1:nrow(ranzig.probes)){
  rowInKrank <- which(ranzig.krank[, "OHR"] == ranzig.probes[x, "ID"])
  ranzig.probes[x, "MastitisInHerde"]  <- length(rowInKrank)
  if(length(rowInKrank) > 0){
    ranzig.probes[x, "Lactation"] <- paste0(unique(ranzig.krank[rowInKrank, "LAKTATION"]), collapse=",")
  }
}
#table(ranzig.probes[,"MastitisInHerde"])
#ranzig.probes[grep("1", ranzig.probes[,"Lactation"]),]

# Simulation
ranzig.probes <- cbind(ranzig.probes, SimMarker = NA)
set.seed(30)
ranzig.probes[,"SimMarker"] <- sample(c("AA", "AT", "TT"), nrow(ranzig.probes), replace=TRUE)

sick <- which(ranzig.probes[,"MastitisInHerde"] > 0)
healthy <- which(ranzig.probes[,"MastitisInHerde"] == 0)

contTable <- rbind(sick = table(ranzig.probes[sick, "SimMarker"]), healthy = table(ranzig.probes[healthy, "SimMarker"]))

anova(lm(ranzig.probes[, "MastitisInHerde"] ~ ranzig.probes[, "SimMarker"]))
