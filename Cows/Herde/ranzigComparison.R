# Herde matching script Junling
#
library(foreign)
setwd("D:/Edrive/Cow/HERDE/NeuDaten08/")

D2018_2  <- cbind(F = "2018_2", read.dbf("Herde/2018_02/BESTAND.DBF"))
D2018_7  <- cbind(F = "2018_7", read.dbf("Herde/2018_07/BESTAND.DBF"))
D2018_7b <- cbind(F = "2018_7b", read.dbf("Herde/2018_07b/BESTAND.DBF"))
D2018_L  <- cbind(F = "2018_L", read.dbf("Herde/HW_91237082_180830/BESTAND.DBF"))