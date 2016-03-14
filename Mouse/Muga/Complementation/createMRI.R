
createMRItable <- function(MRIdata, description){
  dates <- strsplit(unlist(lapply(strsplit(as.character(MRIdata[,"TimeDateDura"]),";"),"[",1))," ")                       # Get the malformed dates

  MRIdata[,"TimeDateDura"] <- unlist(lapply(dates, function(x){                                                           # Transform to DD/MM/YY
    monthNumber <- which(months[,1]==x[2])
    paste(gsub(",","",x[3]), monthNumber, x[4],sep="/")
  }))

  for(x in 1:nrow(MRIdata)){
    bDay <- as.character(description[as.character(MRIdata[x,"Label"]), "W-dat"])                                                          # Date of birth
    mDay <- as.character(MRIdata[x, "TimeDateDura"])                                                                                    # Date of measurement
    daysDiff <- as.numeric(round(difftime(strptime(mDay, format = "%d/%m/%Y"), strptime(bDay, format = "%d.%m.%Y"), units="days")))
    if(length(daysDiff) == 0) daysDiff<- 666                                                                                            # If one of the dates is missing use 666
    cat(bDay, mDay, daysDiff, "\n")
    MRIdata[x, "Age"] <- daysDiff
  }

  animals <- unique(as.character(MRIdata[,"Label"]))
  timepoints <- unique(as.character(MRIdata[,"Age"]))

  fat <- matrix(NA, length(animals), length(timepoints), dimnames=list(animals, timepoints))
  lean <- matrix(NA, length(animals), length(timepoints), dimnames=list(animals, timepoints))

  for(tp in timepoints){
    for(animal in animals){
      ii <- which(MRIdata[,"Label"] == animal & MRIdata[,"Age"] == tp)
      if(length(ii) > 0){
        fat[animal, tp] <- mean(MRIdata[ii,"Fat"])
        lean[animal, tp] <- mean(MRIdata[ii,"Lean"])
      }
    }
  }
  return(list(fat, lean))
}