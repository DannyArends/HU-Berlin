setwd("D:/Edrive/Horse/DNA/SNPchip_Monika_Feb19")
mdata <- readLines("HorseBreeds_SNPchip.txt")

splitted <- lapply(mdata, strsplit, "\t")

header <- splitted[1:3]
splitted <- splitted[-c(1:3)]

nRow <- length(splitted)
nCol <- length(splitted[[1]][[1]]) - 5
mData <- matrix(NA, nRow, nCol)
nItems <- length(splitted[[1]][[1]])

for(x in 1:nRow){
  mData[x,] <- splitted[[x]][[1]][6:nItems]
}

horses <- unique(header[[2]][[1]][6:nItems])
myData <- matrix(NA, nRow, length(horses))
colnames(myData) <- horses
rownames(myData) <- unlist(lapply(lapply(splitted, unlist), "[",1))

colnames(mData) <- header[[2]][[1]][6:nItems]
for(cX in seq(1, 94, 2)){
  myData[,colnames(mData)[cX]] <- paste0(mData[,cX],mData[,(cX+1)])
}
for(cX in 95:ncol(mData)){
  myData[,colnames(mData)[cX]] <- mData[,cX]
}

myMap <- cbind(unlist(lapply(lapply(splitted, unlist), "[",1)),
               unlist(lapply(lapply(splitted, unlist), "[",2)),
               unlist(lapply(lapply(splitted, unlist), "[",3)),
               unlist(lapply(lapply(splitted, unlist), "[",4)),
               unlist(lapply(lapply(splitted, unlist), "[",5)))

colnames(myMap) <- c("ID", "Index", "Chr", "Position", "SNP")
write.table(myData, "genotypes.txt", sep = "\t", quote=FALSE)
write.table(myMap, "map.txt", sep = "\t", quote=FALSE, row.names=FALSE)

genotypes <- read.csv("genotypes.txt", sep = "\t", row.names=1, na.strings=c("", "NA", "--"), colClasses="character")
map <-  read.csv("map.txt", sep = "\t", row.names = 1)

noSeg <- c()
numgeno <- matrix(NA, nrow(genotypes), ncol(genotypes), dimnames=list(rownames(genotypes), colnames(genotypes)))
for(x in 1:nrow(genotypes)) {
  alleles <- as.character(na.omit(unique(unlist(strsplit(as.character(genotypes[x,]), "")))))
  if(length(alleles) == 2) {
    HomA <- paste0(alleles[1],alleles[1])
    Het0 <- paste0(alleles[1],alleles[2])
    Het1 <- paste0(alleles[2],alleles[1])
    HomB <- paste0(alleles[2],alleles[2])
    numgeno[x, genotypes[x,] == HomA] <- 0
    numgeno[x, genotypes[x,] == Het0] <- 1
    numgeno[x, genotypes[x,] == Het1] <- 1
    numgeno[x, genotypes[x,] == HomB] <- 2
  }else{
    noSeg <- c(noSeg, x)
  }
}
numgeno <- numgeno[-noSeg, ]
write.table(numgeno, "numgeno.txt", sep = "\t", quote=FALSE)
map <- map[-noSeg, ]
write.table(map, "numgeno_map.txt", sep = "\t", quote=FALSE)

plot(hclust(dist(t(numgeno), method="manhattan")),hang=-1)
