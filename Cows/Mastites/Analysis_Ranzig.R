setwd("D:/Edrive/Cow/HERDE")
mydata <- read.table("data.txt", sep="\t", header=TRUE)
isMastites = which(grepl("1.13", mydata[,"STAUFEN"], fixed=TRUE))
mastitesData = mydata[isMastites, ]

# Covert dates to R format, and store as new column
RdateObjects = as.Date(as.character(mastitesData[, "DATUM"]), format="%m/%d/%Y")
mastitesData = cbind(mastitesData, "RDate" = RdateObjects)

## Order by date
ordering = order(mastitesData[, "RDate"])
mastitesData = mastitesData[ordering, ]

mastitesData <- cbind(mastitesData, year = unlist(lapply(strsplit(as.character(mastitesData[,"RDate"]), "-"), "[", 1)))

# Select all the entries not between 2009 and 2017
outOfTimeRange = which(!(as.numeric(as.character(mastitesData[, "year"])) > 2009 
                         & as.numeric(as.character(mastitesData[, "year"])) < 2017))
mastitesData = mastitesData[-outOfTimeRange, ]

## Here we can subset our data to do this per year
in2010 <- which(grepl("2010", mastitesData[,"RDate"]))
in2011 <- which(grepl("2011", mastitesData[,"RDate"]))
in2012 <- which(grepl("2012", mastitesData[,"RDate"]))
in2013 <- which(grepl("2013", mastitesData[,"RDate"]))
in2014 <- which(grepl("2014", mastitesData[,"RDate"]))
in2015 <- which(grepl("2015", mastitesData[,"RDate"]))
in2016 <- which(grepl("2016", mastitesData[,"RDate"]))

generateInfo = function(mastitesData){
  mastitesTable = table(as.character(mastitesData[,"OHR"]))
  mastitesInfo <- NULL
  mastitesEntry <- vector("list", length(names(mastitesTable)))
  names(mastitesEntry) = names(mastitesTable)
  for(cow in names(mastitesTable)){
    entriesForCow = which(mastitesData[, "OHR"] == cow)
    #nEntriesCow = length(mastitesData[entriesForCow, "RDate"])
    mastitesEntry[[cow]] = mastitesData[entriesForCow, "RDate"]
    nDiagnosis = length(mastitesEntry[[cow]])
    if(nDiagnosis > 1){
      diagnosis = c("N", rep(NA, nDiagnosis-1))
      diagnosisDate = mastitesEntry[[cow]][1]
      for(x in 2:nDiagnosis){
        if((mastitesEntry[[cow]][x] - diagnosisDate) < 14){
          diagnosis[x] <- "F"
        }else{
          diagnosis[x] <- "N"
        }
        diagnosisDate <- mastitesEntry[[cow]][x]
      }
      mastitesEntry[[cow]] <- rbind(as.character(mastitesEntry[[cow]]), diagnosis)
      nMastites = length(which(diagnosis == "N"))
      mastitesInfo = rbind(mastitesInfo, c(cow, nDiagnosis, nMastites))
    }else{
      mastitesInfo = rbind(mastitesInfo, c(cow, 1, 1))
    }
    #cat("cow",cow,"\n")
  }
  colnames(mastitesInfo) <- c("OHR", "nDiagnosis", "nMastites")
  return(mastitesInfo)
}

mastitesInfo = generateInfo(mastitesData)
sum(as.numeric(mastitesInfo[,"nDiagnosis"])) / nrow(mastitesInfo)
sum(as.numeric(mastitesInfo[,"nMastites"])) / nrow(mastitesInfo)

mastitesInfo = generateInfo(mastitesData[c(in2010,in2011, in2012, in2013, in2014, in2015, in2016), ])
sum(as.numeric(mastitesInfo[,"nDiagnosis"])) / nrow(mastitesInfo)
sum(as.numeric(mastitesInfo[,"nMastites"])) / nrow(mastitesInfo)

year = 2010
for(individuals in list(in2010, in2011, in2012, in2013, in2014, in2015, in2016)){
  mastitesInfo = generateInfo(mastitesData[individuals, ])
  cat("Year: ", year)
  cat(", Cows: ", nrow(mastitesInfo))
  cat(", Diagnosis: ", sum(as.numeric(mastitesInfo[,"nDiagnosis"])))
  cat(", Mastites: ", sum(as.numeric(mastitesInfo[,"nMastites"])), "\n")
  year = year + 1
}

## Per lactation
inL1 <- which(mastitesData[,"LAKTATION"] == 1)
inL2 <- which(mastitesData[,"LAKTATION"] == 2)
inL3 <- which(mastitesData[,"LAKTATION"] == 3)
inL4 <- which(mastitesData[,"LAKTATION"] == 4)
inL5plus <- which(mastitesData[,"LAKTATION"] > 4)

lactation = 1
for(individuals in list(inL1, inL2, inL3, inL4, inL5plus)){
  mastitesInfo = generateInfo(mastitesData[individuals, ])
  cat("Lactation: ", lactation)
  cat(", Cows: ", nrow(mastitesInfo))
  cat(", Diagnosis: ", sum(as.numeric(mastitesInfo[,"nDiagnosis"])))
  cat(", Mastites: ", sum(as.numeric(mastitesInfo[,"nMastites"])), "\n")
  lactation = lactation + 1
}

#Load the Lactation / Kalbung data
lactData= read.table("data_lactation.txt", sep="\t", header=TRUE)
lactData[1:10,]

# Covert dates to R format, and store as new column
KdateObjects = as.Date(as.character(lactData[, "KALBUNG"]), format="%m/%d/%Y")
lactData = cbind(lactData, "KDate" = KdateObjects)

#Order by date of kalving
## Order by date
ordering = order(lactData[, "KDate"])
lactData = lactData[ordering, ]

### Analuze the first lactation for when a novel mastites incidence occurs
mastDataL1 = mastitesData[inL1, ]
cowsInLactation1 = unique(as.character(mastDataL1[, "OHR"]))
lactDataL1 = lactData[which(lactData[, "OHR"] %in% cowsInLactation1),]

allMastitesDates <- NULL
for(cow in cowsInLactation1){
  cat(cow)
  lactStartDate = lactDataL1[which(lactDataL1[,"OHR"] == cow), "KDate"]
  cat(" ", as.character(lactStartDate), " ")
  diagDates = mastDataL1[which(mastDataL1[,"OHR"] == cow), "RDate"]
  cat(diagDates - lactStartDate, " ")    
  dateDiffs = diagDates - lactStartDate
  datesNewMastites <- NULL
  if(length(dateDiffs) > 1){
    diagnosisDate <- dateDiffs[1]
    datesNewMastites <- c(datesNewMastites, diagnosisDate)
    for(x in 2:length(dateDiffs)){
      if((dateDiffs[x] - diagnosisDate) >= 14){
        datesNewMastites <- c(datesNewMastites, dateDiffs[x])
      }
      diagnosisDate <- dateDiffs[x]
    }
  }else{
    datesNewMastites <- c(datesNewMastites, dateDiffs[1])
  }
  cat(" - ", datesNewMastites, " ")
  allMastitesDates <- c(allMastitesDates, datesNewMastites)
  cat("\n")
}

histData = hist(allMastitesDates, breaks = seq(0, 700, 20))
percentageData = 100 * histData$counts / length(cowsInLactation1)

plot(c(0, 700), c(0, 60), t = 'n', las=2, ylab="Percentage", xlab="Days after calfing")
rect(seq(0, 693, 20), rep(0, length(percentageData)), seq(20, 700, 20), percentageData)
