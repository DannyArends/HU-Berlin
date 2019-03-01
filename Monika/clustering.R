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
rownames(myData) <- gsub("_", "-", unlist(lapply(lapply(splitted, unlist), "[",1 )))

colnames(mData) <- header[[2]][[1]][6:nItems]
for(cX in seq(1, 94, 2)){
  myData[,colnames(mData)[cX]] <- paste0(mData[,cX],mData[,(cX+1)])
}
for(cX in 95:ncol(mData)){
  myData[,colnames(mData)[cX]] <- mData[,cX]
}

myMap <- cbind(gsub("_", "-", unlist(lapply(lapply(splitted, unlist), "[",1 ))),
               unlist(lapply(lapply(splitted, unlist), "[",2)),
               unlist(lapply(lapply(splitted, unlist), "[",3)),
               unlist(lapply(lapply(splitted, unlist), "[",4)),
               unlist(lapply(lapply(splitted, unlist), "[",5)))

colnames(myMap) <- c("ID", "Index", "Chr", "Position", "SNP")
write.table(myData, "genotypes.txt", sep = "\t", quote=FALSE)
write.table(myMap, "map.txt", sep = "\t", quote=FALSE, row.names=FALSE)

genotypes <- read.csv("genotypes.txt", sep = "\t", row.names=1, na.strings=c("", "NA", "--"), colClasses="character")
map <-  read.csv("map.txt", sep = "\t", row.names = 1)

### Add peterson data
setwd("D:/Edrive/Horse/DNA/Petersen2013/raw")
for(mfile in list.files(path = ".", pattern = "csv.gz")){
  splitted <- strsplit(readLines(gzfile(mfile)), ",")
  markers <- unlist(lapply(splitted, "[", 1))
  horses <- unlist(lapply(splitted, "[", 2))
  GTs <- paste0(unlist(lapply(splitted, "[", 5)), unlist(lapply(splitted, "[", 6)))
  for (horse in na.omit(unique(horses))) {
    genotypes <- cbind(genotypes, cHorse = NA)
    colnames(genotypes)[which(colnames(genotypes) == "cHorse")] <- horse
    cHorseData <- which(horses == horse)
    cHorseGTs <- GTs[cHorseData]
    names(cHorseGTs) <- markers[cHorseData]
    genotypes[, horse] <- cHorseGTs[rownames(genotypes)]
    cat("Done horse: ", horse, "\n")
  }
}
genotypes[genotypes == "--"] <- NA
setwd("D:/Edrive/Horse/DNA/SNPchip_Monika_Feb19")
write.table(genotypes, "genotypes_peterson_merged.txt", sep = "\t", quote=FALSE)

opposite <- function(x){
  if(any(is.na(x))) return("")
  ret <- NULL
  for(e in x){
    if(e == "A") ret <- c(ret, "T")
    if(e == "C") ret <- c(ret, "G")
    if(e == "G") ret <- c(ret, "C")
    if(e == "T") ret <- c(ret, "A")
  }
  return(ret)
}

# Fix direction of alleles in Peterson data
unfixable <- c()
for(x in 1:nrow(genotypes)) {
  alleles_monika <- sort(as.character(na.omit(unique(unlist(strsplit(as.character(genotypes[x,1:92]), ""))))))
  alleles_peterson <- sort(as.character(na.omit(unique(unlist(strsplit(as.character(genotypes[x,93:ncol(genotypes)]), ""))))))
  if(length(alleles_monika) == 2 && length(alleles_peterson) == 2 && !all(alleles_monika == alleles_peterson)){
    nAlleles <- opposite(alleles_peterson)
    genotypes[x, 93:ncol(genotypes)] <- gsub(alleles_peterson[2], nAlleles[2], gsub(alleles_peterson[1], nAlleles[1], genotypes[x,93:ncol(genotypes)]))
  }else{
    cat(x, " ", length(alleles_monika), " ", length(alleles_peterson), "\n")
  }
}
write.table(genotypes, "genotypes_peterson_merged_and_flipped.txt", sep = "\t", quote=FALSE)


noSeg <- c()
numgeno <- matrix(NA, nrow(genotypes), ncol(genotypes), dimnames=list(rownames(genotypes), colnames(genotypes)))
for(x in 1:nrow(genotypes)) {
  alleles <- sort(as.character(na.omit(unique(unlist(strsplit(as.character(genotypes[x,]), ""))))))
  if(length(alleles) == 2) {
    HomA <- paste0(alleles[1],alleles[1])
    Het0 <- paste0(alleles[1],alleles[2])
    Het1 <- paste0(alleles[2],alleles[1])
    HomB <- paste0(alleles[2],alleles[2])
    numgeno[x, genotypes[x,] == HomA] <- 0
    numgeno[x, genotypes[x,] == Het0] <- 1
    numgeno[x, genotypes[x,] == Het1] <- 1
    numgeno[x, genotypes[x,] == HomB] <- 2
  } else if(length(alleles) < 2) {
    cat("No seg at",x,"\n")
    noSeg <- c(noSeg, x)
  } else {
    noSeg <- c(noSeg, x)
    cat("[WARNING] ", x,"more then 2 alleles present (", length(alleles), ") shouldnt have happened after flipping:",alleles,"\n")
  }
}


numgeno <- numgeno[-noSeg, ]
write.table(numgeno, "numgeno_peterson_merged.txt", sep = "\t", quote=FALSE)
map <- map[-noSeg, ]
write.table(map, "numgeno_peterson_merged_map.txt", sep = "\t", quote=FALSE)

setwd("D:/Edrive/Horse/DNA/SNPchip_Monika_Feb19")
numgeno <- read.table("numgeno_peterson_merged.txt", sep = "\t")
map <- read.table("numgeno_peterson_merged_map.txt", sep = "\t")

monika <- 1:92
caspian <- grep("CS_", colnames(numgeno)) # Caspian
colnames(numgeno)[caspian] <- paste0("Caspian ", 1:length(caspian))
tuva <- grep("Tu", colnames(numgeno)) # Tuva
colnames(numgeno)[tuva] <- paste0("Tuva ", 1:length(tuva))
mongolian <- grep("_MON", colnames(numgeno)) # Mongolian
colnames(numgeno)[mongolian] <- paste0("Mongolian ", 1:length(mongolian))
exmoor <- grep("EX_", colnames(numgeno)) # Exmoor
colnames(numgeno)[exmoor] <- paste0("Exmoor ", 1:length(exmoor))
thoroughbred <- grep("_TB", colnames(numgeno)) # Thoroughbred
colnames(numgeno)[thoroughbred] <- paste0("Thoroughbred ", 1:length(thoroughbred))
arabian <- grep("ARR", colnames(numgeno)) # Arabian
colnames(numgeno)[arabian] <- paste0("Arabian ", 1:length(arabian))
akhalteke <- grep("AH_", colnames(numgeno)) # Akhal_Teke
colnames(numgeno)[akhalteke] <- paste0("AkhalTeke ", 1:length(akhalteke))

# Columns of interest
cOI <- c(monika, caspian, tuva, mongolian, exmoor, thoroughbred, arabian, akhalteke)

n <- 2500
r <- 50
sS <- (nrow(numgeno) / n)

distances <- vector("list", r)
for(x in 1:r){
  distances[[x]] <- dist(t(numgeno[sample(nrow(numgeno), n), cOI]), method="manhattan")
  cat("Done", x, "/", r, "\n")
}

dSum <- (distances[[1]] * sS)
for(x in 2:r){
  dSum <- dSum + (distances[[x]] *sS)
}
dAvg <- dSum / r

library(dendextend) # library(circlize)
library(RColorBrewer)

clustering <- hclust(dAvg)
ordering <- clustering$labels[clustering$order]
mdendrogram <- as.dendrogram(clustering)

mcolors = brewer.pal(7, "Paired")
names(mcolors) <- c("Caspian", "Tuva", "Mongolian", "Exmoor", "Thoroughbred", "Arabian", "AkhalTeke")

labelCol <- function(x){
  if(is.leaf(x)){
    label <- strsplit(attr(x, "label"), " ")[[1]][1]
    attr(x, "nodePar") <- list(lab.col = mcolors[label], pch=NA)  # Set the label color based on the strain
    #attr(x, "label") <- label
  }
  return(x)
}

dendrocol <- dendrapply(mdendrogram, labelCol)

op <- par(cex=0.7)
plot(dendrocol)
circlize_dendrogram(dendrocol, labels_track_height=NA)
plot(dendrocol)