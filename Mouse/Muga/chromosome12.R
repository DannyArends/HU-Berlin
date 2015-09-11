# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sept, 2015
# first written Sept, 2015

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

setwd("E:/Mouse/DNA/MegaMuga/")
mdata <- read.table("data_chromosome12_complete.txt", sep = "\t",header=TRUE, stringsAsFactor=FALSE,colClasses=c("character"))

mdata[which(!as.numeric(mdata[,"TierID"]) < 4000),"TierID"] <- paste0("333", mdata[which(!as.numeric(mdata[,"TierID"]) < 4000),"TierID"]) # Children add 333
mdata[which(as.numeric(mdata[,"TierID"]) < 4000),"TierID"] <- paste0("666", mdata[which(as.numeric(mdata[,"TierID"]) < 4000),"TierID"]) # Children add 666

#unlist(lapply(strsplit(as.character(mdata[,"TierID"]),""),length))  # Length of 5 or tierID > 4000 -> Add 333
phenotypedata <- phenotypedata[match(mdata[,"TierID"], rownames(phenotypedata)),]
 
combined <- cbind(mdata, phenotypedata)
individuals <- rownames(combined)

markerIDs <- c("KM.030", "KM.044", "KM.045", "KM.046")

for(x in individuals){
  vID <- as.character(combined[x,"Vater"])
  mID <- as.character(combined[x,"Mutter"])
  if(!is.na(vID) && !is.na(mID)){
    datamatrix <- NULL
    for(m in markerIDs){
      vGeno <- strsplit(combined[vID, m],"")[[1]]
      if(length(vGeno) != 2) vGeno <- c("?","?")
      mGeno <- strsplit(combined[mID, m],"")[[1]]
      if(length(mGeno) != 2) mGeno <- c("?","?")
      cGeno <- strsplit(combined[x, m],"")[[1]]
      if(length(cGeno) != 2) cGeno <- c("?","?")
    
      datamatrix <- rbind(datamatrix, c("M", m, vGeno, mGeno, cGeno))
    }
    input  <- paste0("Analysis/rephase/", x, ".bgl")
    output <- paste0("Analysis/rephase/", x)

    write.table(datamatrix, file = input, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)
    system(paste0("java -Xmx1000m -jar beagle.jar redundant=true trios=", input, " missing=? out=", output))
  }
}

phaseddata <- matrix(NA, length(individuals),length(markerIDs), dimnames=list(individuals, markerIDs))

for(x in individuals){
  fname <- paste0("Analysis/rephase/", x,".",x,".bgl.phased.gz")
  if(file.exists(fname)){
    pdata <- read.table(gzfile(paste0("Analysis/rephase/", x,".",x,".bgl.phased.gz")), header=TRUE, colClasses="character", row.names=2)
    phaseddata[x, rownames(pdata)] <-  apply(pdata[,6:7],1,function(x){paste0(x[1],x[2]); })
  }else{
    cat("individual",x,"no data\n")
  }
}



