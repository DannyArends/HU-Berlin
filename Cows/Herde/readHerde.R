
library(foreign)
setwd("D:/Edrive/Cow/HERDE/NeuDaten08/")

path <- "Herde/HW_91237082_180830/"

dbffiles <- list.files(path, "*.DBF")
herdedata <- vector("list", length(dbffiles))
names(herdedata) <- tolower(dbffiles)
for(dbffile in dbffiles){
  cat("Loading: ", dbffile, "\n")
  herdedata[tolower(dbffile)]  <- read.dbf(paste0(path, dbffile), as.is = TRUE)
}


