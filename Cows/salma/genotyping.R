setwd("D:/Edrive/Cow/HERDE/Salma")

barcodes <- read.table("barcodes.txt", sep = "\t", header=TRUE)
plate_384 <- read.table("plate_384.txt",sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
plate_96 <- read.table("plate_96.txt",sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
genotypes_384 <- read.table("genotypes_384.txt", sep = "\t", header = TRUE)
genotypes_96 <- read.table("genotypes_96.txt", sep = "\t", header = TRUE)
phenotypes <- read.table("phenotypes.txt", sep="\t", header=TRUE)

genotypes_384 <- genotypes_384[,c("Well.Position", "Call")]
genotypes_384[,"Call"] <- gsub("Heterozygous ", "", gsub("Homozygous ", "", genotypes_384[,"Call"]))
genotypes_384 <- cbind(Sample = NA, genotypes_384)

genotypes_96 <- genotypes_96[,c("Well", "Call")]
genotypes_96[,"Call"] <- gsub("Heterozygous ", "", gsub("Homozygous ", "", genotypes_96[,"Call"]))
genotypes_96 <- cbind(Sample = NA, genotypes_96)
genotypes_96[which(genotypes_96[,"Call"] == ""),] <- NA

bcodeIDS <- gsub("DE", "DE00", as.character(barcodes[,1]))
hasPheno <- bcodeIDS[which(bcodeIDS %in% unique(phenotypes[,"id"]))]

for(x in rownames(plate_384)){
  for(y in colnames(plate_384)){
    wellID <- paste0(x, y)
    gtID <- which(genotypes_384[, "Well.Position"] == wellID)
    genotypes_384[gtID, "Sample"] <- as.character(plate_384[x,y])
  }
}

for(x in rownames(plate_96)){
  for(y in colnames(plate_96)){
    wellID <- paste0(x, y)
    gtID <- which(genotypes_96[, "Well"] == wellID)
    genotypes_96[gtID, "Sample"] <- as.character(plate_96[x,y])
  }
}

genotypes_384 <- cbind(Ohr = NA, genotypes_384)
for(x in 1:nrow(genotypes_384)){
  bCodeID <- which(barcodes[,"Lab.Nr"] == paste0("E", genotypes_384[x, "Sample"]))
  if(length(bCodeID) > 0){
    genotypes_384[x, "Ohr"] <- as.character(barcodes[bCodeID, "Barcode.Ohmark"])
  }else{
    genotypes_384[x, "Ohr"] <- "CTRL"
  }
}
genotypes_384[, "Ohr"] <- gsub("DE", "DE00", genotypes_384[, "Ohr"] )

genotypes_96 <- cbind(Ohr = NA, genotypes_96)
for(x in 1:nrow(genotypes_96)){
  bCodeID <- which(barcodes[,"Lab.Nr"] == paste0("E", genotypes_96[x, "Sample"]))
  if(length(bCodeID) > 0){
    genotypes_96[x, "Ohr"] <- as.character(barcodes[bCodeID, "Barcode.Ohmark"])
  }else{
    genotypes_96[x, "Ohr"] <- "CTRL"
  }
}
genotypes_96[, "Ohr"] <- gsub("DE", "DE00", genotypes_96[, "Ohr"] )

phenotypes <- cbind(phenotypes, GT = NA, Plate = NA)
phenotypes <- cbind(phenotypes, year = unlist(lapply(strsplit(as.character(phenotypes[,"birth.date"]), "/"), "[",3)))

for(x in 1:nrow(phenotypes)){
  gtID <- which(genotypes_384[, "Ohr"] == phenotypes[x, "id"])
  if(length(gtID) > 0){
    phenotypes[x , "GT"] <- unique(genotypes_384[gtID, "Call"])
    phenotypes[x , "Plate"] <- "384"
  }
  gtID <- which(genotypes_96[, "Ohr"] == phenotypes[x, "id"])
  if(length(gtID) > 0){
    phenotypes[x , "GT"] <- unique(genotypes_96[gtID, "Call"])
    phenotypes[x , "Plate"] <- "96"
  }
}

phenotypes <- phenotypes[which(!is.na(phenotypes[,"GT"])),]
phenotypes <- phenotypes[which(phenotypes[,"lactation"] == 2),]

library(heterozygous)
HWE(gsub("/", "",(t(phenotypes[,"GT"]))))

tbl <- table(phenotypes[, c("Erkrankung")])
names(tbl)[1] <- "Keine Mastitis"
barplot(tbl, col=c("darkblue", "darkred"), main="Mastitisverteilung")
text(0.70, 10, paste0("n = ", tbl[1]), col="white")
text(1.90, 10, paste0("n = ", tbl[2]), col="white")

tbl <- table(phenotypes[,"GT"])
barplot(tbl, col=c("red", "green", "blue"), main="Genotypverteilung")
text(0.70,10, paste0("n = ", tbl[1]))
text(1.90,10, paste0("n = ", tbl[2]))
text(3.10,10, paste0("n = ", tbl[3]), col="white")

tbl <- table(phenotypes[, c("Erkrankung", "GT")])
barplot(tbl, col=c("darkblue", "darkred"), main="Mastitis pro Genotyp", ylim=c(0,300))
legend("topleft", c("Keine Mastitis", "Mastitis"), fill=c("darkblue", "darkred"))
text(0.73,10, paste0(round(tbl[2,1] / sum(tbl[,1]) * 100,1), "%"))
text(1.92,97, paste0(round(tbl[2,2] / sum(tbl[,2]) * 100,1), "%"))
text(3.10,290, paste0(round(tbl[2,3] / sum(tbl[,3]) * 100,1), "%"))

chisq.test(tbl)

mtbl <- cbind(tbl[,1] + tbl[,2], tbl[,3])
chisq.test(mtbl)


nAA <- tbl[1,1] * 2 + tbl[1,2]; nGG <- tbl[1,3] * 2 + tbl[1,2]
noM <- nAA / (nAA+nGG)

mtab <- c(nAA, nGG)

nAA <- tbl[2,1] * 2 + tbl[2,2]; nGG <- tbl[2,3] * 2 + tbl[2,2]
yesM <- nAA / (nAA+nGG)

mtab <- rbind(mtab, c(nAA, nGG))
rownames(mtab) <- c("No", "Yes")
colnames(mtab) <- c("A", "G")
chisq.test(mtab)


lm(phenotypes[, "Mkg"] ~ phenotypes[,"GT"])

toGenotype <- hasPheno[!hasPheno %in%  phenotypes[,"id"]]
toGenotype <- gsub("DE00", "DE", toGenotype)

barcodes[which(barcodes[,1] %in% toGenotype),]
write.table(barcodes[which(barcodes[,1] %in% toGenotype),], file= "todo.txt",sep="\t")

