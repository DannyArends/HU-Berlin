setwd("D:/Edrive/Cow/Salma Wachow")
phenotypes <- read.table("phenotypes.txt", sep="\t", header=TRUE)
phenotypes <- cbind(phenotypes, year = unlist(lapply(strsplit(as.character(phenotypes[,"birth.date"]), "/"), "[",3)))

barcodes <- read.table("barcodes.txt", sep = "\t", header=TRUE)
barcodes[,1] <- gsub("DE", "DE00", as.character(barcodes[,1]))
rownames(barcodes) <- barcodes[,"Lab.Nr"]

## SNP 1: rs41257360
rs41257360 <- read.table("Salma_E/R-HAS-M1-5-_rs41257360/rs41257360.txt",sep="\t", header=TRUE)
rs41257360.samples <- read.table("Salma_E/R-HAS-M1-5-_rs41257360/samples.txt",sep="\t", header=TRUE, row.names=1,check.names=FALSE)

rs41257360.data <- cbind(as.character(rs41257360[,"Well.Position"]), unlist(lapply(strsplit(as.character(rs41257360[,"Call"]), " "),"[",2)))
colnames(rs41257360.data) <- c("Well.Position", "GT")

rs41257360.data <- cbind(rs41257360.data, rowID = unlist(lapply(strsplit(rs41257360.data[,1], ""),"[",1)))
rs41257360.data <- cbind(rs41257360.data, colID = as.numeric(substring(rs41257360.data[,1], 2)))

labnrs <- c()
for(x in 1:nrow(rs41257360.data)){
  labnrs <- c(labnrs, as.character(rs41257360.samples[rs41257360.data[x, "rowID"], rs41257360.data[x, "colID"]]))
}

rs41257360.data <- cbind(rs41257360.data, "Lab.Nr" = labnrs, "Barcode.Ohmark" = barcodes[labnrs, "Barcode.Ohmark"])
rs41257360.data[1:10,]

## SNP 2: rs41634110
rs41634110 <- read.table("Salma_E/R-HAS-M3-6-_rs41634110/rs41634110.txt",sep="\t", header=TRUE)
rs41634110.samples <- read.table("Salma_E/R-HAS-M3-6-_rs41634110/samples.txt",sep="\t", header=TRUE, row.names=1,check.names=FALSE)

rs41634110.data <- cbind(as.character(rs41634110[,"Well.Position"]), unlist(lapply(strsplit(as.character(rs41634110[,"Call"]), " "),"[",2)))
colnames(rs41634110.data) <- c("Well.Position", "GT")

rs41634110.data <- cbind(rs41634110.data, rowID = unlist(lapply(strsplit(rs41634110.data[,1], ""),"[",1)))
rs41634110.data <- cbind(rs41634110.data, colID = as.numeric(substring(rs41634110.data[,1], 2)))

labnrs <- c()
for(x in 1:nrow(rs41634110.data)){
  labnrs <- c(labnrs, as.character(rs41634110.samples[rs41634110.data[x, "rowID"], rs41634110.data[x, "colID"]]))
}

rs41634110.data <- cbind(rs41634110.data, "Lab.Nr" = labnrs, "Barcode.Ohmark" = barcodes[labnrs, "Barcode.Ohmark"])
rs41634110.data[1:10,]

## SNP 3: rs29020544
rs29020544 <- read.table("Salma_E/R-HAS-M7-9-_rs29020544/rs29020544.txt",sep="\t", header=TRUE)
rs29020544.samples <- read.table("Salma_E/R-HAS-M7-9-_rs29020544/samples.txt",sep="\t", header=TRUE, row.names=1,check.names=FALSE)

rs29020544.data <- cbind(as.character(rs29020544[,"Well.Position"]), unlist(lapply(strsplit(as.character(rs29020544[,"Call"]), " "),"[",2)))
colnames(rs29020544.data) <- c("Well.Position", "GT")

rs29020544.data <- cbind(rs29020544.data, rowID = unlist(lapply(strsplit(rs29020544.data[,1], ""),"[",1)))
rs29020544.data <- cbind(rs29020544.data, colID = as.numeric(substring(rs29020544.data[,1], 2)))

labnrs <- c()
for(x in 1:nrow(rs29020544.data)){
  labnrs <- c(labnrs, as.character(rs29020544.samples[rs29020544.data[x, "rowID"], rs29020544.data[x, "colID"]]))
}

rs29020544.data <- cbind(rs29020544.data, "Lab.Nr" = labnrs, "Barcode.Ohmark" = barcodes[labnrs, "Barcode.Ohmark"])
rs29020544.data[1:10,]

## SNP 4: rs41629005
rs41629005 <- read.table("Salma_E/R-HAS-M10-8-_rs41629005/rs41629005.txt",sep="\t", header=TRUE)
rs41629005.samples <- read.table("Salma_E/R-HAS-M10-8-_rs41629005/samples.txt",sep="\t", header=TRUE, row.names=1,check.names=FALSE)

rs41629005.data <- cbind(as.character(rs41629005[,"Well.Position"]), unlist(lapply(strsplit(as.character(rs41629005[,"Call"]), " "),"[",2)))
colnames(rs41629005.data) <- c("Well.Position", "GT")

rs41629005.data <- cbind(rs41629005.data, rowID = unlist(lapply(strsplit(rs41629005.data[,1], ""),"[",1)))
rs41629005.data <- cbind(rs41629005.data, colID = as.numeric(substring(rs41629005.data[,1], 2)))

labnrs <- c()
for(x in 1:nrow(rs41629005.data)){
  labnrs <- c(labnrs, as.character(rs41629005.samples[rs41629005.data[x, "rowID"], rs41629005.data[x, "colID"]]))
}

rs41629005.data <- cbind(rs41629005.data, "Lab.Nr" = labnrs, "Barcode.Ohmark" = barcodes[labnrs, "Barcode.Ohmark"])
rs41629005.data[1:10,]


# Combine all SNPs in a matrix

samples <- na.omit(unique(c(rs41629005.data[, "Barcode.Ohmark"], rs29020544.data[, "Barcode.Ohmark"],rs41634110.data[, "Barcode.Ohmark"],rs41257360.data[, "Barcode.Ohmark"])))

genotypes <- matrix(NA, length(samples), 4, dimnames=list(samples, c("rs41629005","rs29020544","rs41634110","rs41257360")))

for(x in 1:nrow(rs41629005.data)){
  if(!is.na(rs41629005.data[x, "Barcode.Ohmark"])) genotypes[rs41629005.data[x, "Barcode.Ohmark"],"rs41629005"] <- rs41629005.data[x, "GT"]
}
for(x in 1:nrow(rs29020544.data)){
  if(!is.na(rs29020544.data[x, "Barcode.Ohmark"])) genotypes[rs29020544.data[x, "Barcode.Ohmark"],"rs29020544"] <- rs29020544.data[x, "GT"]
}
for(x in 1:nrow(rs41634110.data)){
  if(!is.na(rs41634110.data[x, "Barcode.Ohmark"])) genotypes[rs41634110.data[x, "Barcode.Ohmark"],"rs41634110"] <- rs41634110.data[x, "GT"]
}
for(x in 1:nrow(rs41257360.data)){
  if(!is.na(rs41257360.data[x, "Barcode.Ohmark"])) genotypes[rs41257360.data[x, "Barcode.Ohmark"],"rs41257360"] <- rs41257360.data[x, "GT"]
}


lactation <- phenotypes[which(phenotypes[,"lactation"] == 5),]
rownames(lactation) <- lactation[,"id"]

npheno <- c("AFC", "firstMkg", "Mkg", "Mkg100", "fatkg100", "fatkg", "proteinkg1", "proteinkg", "countMast")
pvalues <- c()
for(pheno in npheno){
  res <- apply(genotypes, 2, function(marker){
    anova(lm(lactation[rownames(genotypes),pheno] ~ marker))[[5]][1]
  })
  pvalues <- rbind(pvalues, res)
}

rownames(pvalues) <- npheno
pvalues
