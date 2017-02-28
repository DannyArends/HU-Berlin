
### Load in the files from Uni Giessen
setwd("D:/Ddrive/GiesenFugato")
files <- c(paste0("776/", dir("776")), paste0("842/", dir("842")))
results <- vector("list", length(files))
for(i in 1:length(files)) {
  f <- files[i]
  dataf <- read.table(f, skip=10, sep="\t", header=TRUE)
  individualID <- as.character(unique(dataf[,2]))
  GTs <- apply(dataf[,c("Allele1...Top","Allele2...Top")], 1, function(x){
    paste0(unlist(as.character(x)),collapse="")
  })
  names(GTs) <- dataf[,"SNP.Name"]
  GTs[which(dataf[,"GT.Score"] <= 0.7)] <- "--"
  results[[i]] <- list(f, individualID, GTs)
}

indnames <- unlist(lapply(results,"[",2))
probenames <- names(results[[1]][[3]])

resultsM <- matrix(NA, length(results[[1]][[3]]), length(indnames), dimnames=list(probenames, indnames))
x <- lapply(results, function(x){
  resultsM[,x[[2]]] <<- x[[3]][probenames]
})

write.table(resultsM, file="GT_UniGiessen.txt", sep="\t", quote=FALSE)

resultsM <- read.table("GT_UniGiessen.txt", header=TRUE)

ann776 <- read.csv("Samples_Table_LOMs_776.csv", header=TRUE, sep = ";",check.names=FALSE, colClasses='character')
ann842 <- read.csv("Samples_Table_LOMs_842.csv", header=TRUE, sep = ";",check.names=FALSE, colClasses='character')

ann_all <- rbind(ann776, ann842[,-c(10,11)])
annotation <- ann_all[, "LOM"]
names(annotation) <- apply(ann_all[, c("Plattencode", "Positionscode")],1,paste0, collapse="_")

colnames(resultsM) <- annotation[colnames(resultsM)]
resultsM[1:10,1:10]
write.table(resultsM, file="GT_UniGiessen_LOM.txt", sep="\t", quote=FALSE)

