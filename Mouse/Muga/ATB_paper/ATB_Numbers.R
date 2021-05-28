# Load the genotype call data
setwd("D:/Edrive/Mouse/DNA/MegaMuga/")
counts <- read.table(paste0("TRDredo/TransmissionBias_annotated_28.txt"))
regionsMat <- read.table("TRDredo/regions_matp0.01.txt", sep="\t", header=TRUE)
for(r in 1:nrow(regionsMat)){
  sM <- as.character(regionsMat[r, "Flanking.Marker"])
  eM <- as.character(regionsMat[r, "Flanking.Marker.1"])
  sI <- which(rownames(counts) == sM)
  eI <- which(rownames(counts) == eM)
  countsInRegion <- counts[sI:eI,]
  
  top <- countsInRegion[which.min(countsInRegion[, "pMat"])[1],]
  
  swap <- which(countsInRegion[, "bfmi"] == "B")
  
  M12 <- countsInRegion[swap, "M12"]
  M21 <- countsInRegion[swap, "M21"]
  countsInRegion[swap, "M21"] <- M12
  countsInRegion[swap, "M12"] <- M21

  means <- round(apply(countsInRegion[, c("M21", "M12")], 2, mean, na.rm=TRUE),1)
  sds <- round(apply(countsInRegion[, c("M21", "M12")], 2, sd, na.rm=TRUE),1)
  cat("Mat_R",r," ",dim(countsInRegion)[1], ",",  means[1], " (", sds[1], "),",  means[2], " (", sds[2], "),", rownames(top), ", ", top[, "Pos"], ", ", top[,"M12"], ", ", top[,"M21"] , "\n", collapse="", sep ="")
}

# For region 22 and 26 fix them by hand
mean(apply(countsInRegion[, c("M12", "M21")],1,max))
sd(apply(countsInRegion[, c("M12", "M21")],1,max))

mean(apply(countsInRegion[, c("M12", "M21")],1,min))
sd(apply(countsInRegion[, c("M12", "M21")],1,min))

regionsPat <- read.table("TRDredo/regions_patp0.01.txt", sep="\t", header=TRUE)
for(r in 1:nrow(regionsPat)){
  sM <- as.character(regionsPat[r, "Flanking.Marker"])
  eM <- as.character(regionsPat[r, "Flanking.Marker.1"])
  sI <- which(rownames(counts) == sM)
  eI <- which(rownames(counts) == eM)
  countsInRegion <- counts[sI:eI,]
  
  top <- countsInRegion[which.min(countsInRegion[, "pPat"])[1],]

  
  swap <- which(countsInRegion[, "bfmi"] == "B")

  P12 <- countsInRegion[swap, "P12"]
  P21 <- countsInRegion[swap, "P21"]
  countsInRegion[swap, "P21"] <- P12
  countsInRegion[swap, "P12"] <- P21

  means <- round(apply(countsInRegion[, c("P21", "P12")], 2, mean, na.rm=TRUE),1)
  sds <- round(apply(countsInRegion[, c("P21", "P12")], 2, sd, na.rm=TRUE),1)
  cat("Pat_R",r," ",dim(countsInRegion)[1], ",",  means[1], " (", sds[1], "),",  means[2], " (", sds[2], "),", rownames(top), ", ", top[, "Pos"], ", ", top[,"P12"], ", ", top[,"P21"] , "\n", collapse="", sep ="")
}

# For region 15 and 25 fix them by hand
mean(apply(countsInRegion[, c("P12", "P21")],1,max))
sd(apply(countsInRegion[, c("P12", "P21")],1,max))

mean(apply(countsInRegion[, c("P12", "P21")],1,min))
sd(apply(countsInRegion[, c("P12", "P21")],1,min))
