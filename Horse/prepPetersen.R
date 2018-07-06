
setwd("D:/Edrive/Horse/DNA/Petersen2013/raw")

files <- dir()

alldata <- NULL
for(fp in files){
  mdata <- read.csv(gzfile(fp),header=FALSE, colClasses="character")
  breed <- gsub("_", " ", gsub("_Raw.csv.gz", "", fp))
  markers <- unique(mdata[,1])
  samples <- unique(mdata[,2])
  rawdata <- matrix(NA, length(markers), length(samples), dimnames=list(markers, samples))
  for(s in samples){
    msubset <- mdata[which(mdata[,2] == s),c(1,5,6)]
    gts <- paste0(msubset[,2], "/", msubset[,3])
    gts <- gsub("-/-", "./.", gts)
    rawdata[msubset[,1],s] <- gts
  }
  colnames(rawdata) <- paste0(breed, "_", gsub("_", "-", colnames(rawdata)))
  if(is.null(alldata)){
    alldata <- rawdata
  }else{
    alldata <- cbind(alldata, rawdata[rownames(alldata),])
  }
}


breeds <- unlist(lapply(strsplit(colnames(alldata),"_"),"[",1))
hnames <- unlist(lapply(strsplit(colnames(alldata),"_"),"[",2))
indinfo <- cbind(Name = hnames, Breed = breeds)

colnames(alldata) <- hnames

setwd("D:/Edrive/Horse/DNA/Petersen2013/original")

mapdata <- read.csv("mapinfo.txt",sep="\t",row.names=2)
mapdata <- mapdata[rownames(alldata),]
mapdata <- mapdata[-which(mapdata[,"Chr"] == "X"),]
mapdata <- mapdata[sort(mapdata[,"Position"], index.return=TRUE)$ix,]

mapdataordered <- NULL
for(chr in 1:31){
  mapdataordered <- rbind(mapdataordered, mapdata[which(mapdata[,"Chr"] == chr),])
}
mapdata <- mapdataordered

sumpos <- mapdata[,"Position"]
pos.add <- 0
for(chr in 1:31){
  chr.length <- max(mapdata[which(mapdata[,"Chr"]==chr),"Position"],na.rm=TRUE)
  sumpos[which(mapdata[,"Chr"]==chr)] <- sumpos[which(mapdata[,"Chr"]==chr)] + pos.add
  pos.add <- pos.add + chr.length
}

mapdata <- cbind(mapdata, AbsPosition = sumpos)
mapdata <- mapdata[sort(sample(nrow(mapdata), 4000)),]
plot(mapdata[,"AbsPosition"])

alldata <- alldata[rownames(mapdata),]

map <- mapdata[,c(2,3,6)]
colnames(map) <- c("Chromosome", "Position", "SumPosition")

goodbreeds <- names(which(table(indinfo[,"Breed"]) > 20))
indinfo <- indinfo[which(indinfo[,"Breed"] %in% goodbreeds),]

breeds <- indinfo
genotypes <- alldata[,indinfo[,"Name"]]

setwd("D:/Edrive/Horse/DNA/forR")
write.table(indinfo, "individuals.txt",sep="\t", quote=FALSE, row.names=FALSE)
save(breeds, file="breeds.Rdata", compress="xz")
write.table(map, "map.txt",sep="\t", quote=FALSE)
save(map, file="map.Rdata", compress="xz")
write.table(genotypes, "genotypes.txt",sep="\t", quote=FALSE)
save(genotypes, file="genotypes.Rdata", compress="xz")

