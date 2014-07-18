# Created by Danny Arends (c) 2014
# Create mean values for Excel from the output of the MRI machine

##### Change this to the new folder where you save the data #####
setwd("E:/Jan/GMO-Versuche/MRIscript")				
##### Input file from the machine #####
inputMachine <- read.table("AIL-GT 30_MRI_140708.tab", skip=2, sep="\t", row.names=1, colClasses=c("character"))
##### Input file with the mouse birthdays #####
inputDate    <- read.table("GT30-AIL2_Bday.txt", header=TRUE, colClasses=c("character"))

colnames(inputMachine) <- c("ID", "Fat", "Lean", "FreeWater", "TotalWater", "TimeDate", "Unknown")
inputMachine <- inputMachine[-which(inputMachine[,"ID"] == "  "),] 						# Remove the ones without ID, TODO: might be better to fix them

months <- rbind(c("Jan",31),c("Feb",28),c("Mar",31),c("Apr",30),c("May",31),c("Jun",30),c("Jul",31),c("Aug",31),c("Sep",30),c("Oct",31),c("Nov",30),c("Dec",31))

dates <- strsplit(unlist(lapply(strsplit(inputMachine[,"TimeDate"],";"),"[",1))," ")	# Get the malformed dates

inputMachine[,"TimeDate"] <- unlist(lapply(dates, function(x){							# Transform to DD/MM/YY
  monthNumber <- which(months[,1]==x[3])
  paste(gsub(",","",x[4]), monthNumber, x[5],sep="/")
}))

daysDifferences <- NULL																	# Difference in Date, if unknown use 666
for(x in 1:nrow(inputMachine)){
	ID <- inputMachine[x,"ID"]
	bDay <- inputDate[which(inputDate[,"ID"] == gsub(" ","",ID)),"Bday"]
	mDay <- inputMachine[x, "TimeDate"]
	daysDiff <- as.numeric(round(difftime(strptime(mDay, format = "%d/%m/%Y"), strptime(bDay, format = "%d/%m/%y"), units="days")))
	if(length(daysDiff) == 0) daysDiff<- 666
	daysDifferences <- c(daysDifferences, daysDiff)
	cat(x , daysDiff ,"\n")
}

inputMachine <- cbind(inputMachine,Diff = daysDifferences)

columnnames <- c(paste("Fat", unique(daysDifferences),sep="_"), paste("Lean", unique(daysDifferences),sep="_"))
uniqueIDs <- na.omit(unique(inputMachine[,"ID"]))

outPut <- matrix(NA, length(uniqueIDs), length(columnnames))
colnames(outPut) <- columnnames
rownames(outPut) <- uniqueIDs 

for(ID in uniqueIDs){
	IDrows <-  which(inputMachine[,"ID"]==ID)
	differences <- unique(inputMachine[IDrows,"Diff"])
	for(difference in differences){
		toMeanIDs <- IDrows [IDrows %in% which(inputMachine[,"Diff"] == difference)]
		fatMean <- mean(as.numeric(gsub(" ","", inputMachine[toMeanIDs,"Fat"])))
		leanMean <- mean(as.numeric(gsub(" ","", inputMachine[toMeanIDs,"Lean"])))
		outPut[ID,paste("Fat",difference, sep="_")] <- fatMean
		outPut[ID,paste("Lean",difference, sep="_")] <- leanMean
	}
}

write.table(cbind(rownames(outPut),outPut),"output.txt", sep="\t", row.names=FALSE)		# Save for excel, *just drag it there*
