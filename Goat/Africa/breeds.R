setwd("D:/Edrive/Goat/")

#Adaptmap
adaptmap.annot <- read.table("Adaptmap/breeds.txt", sep = "\t", row.names=1, header=TRUE, quote="\"")
#Algeria
algeria.annot <- read.table("Algeria/algeria.annot.txt", sep = "\t", quote="'", header=TRUE, row.names = 1)
#Botswana
ethiopia.annot <- read.table("Ethiopia&Cameron/breeds.txt", sep = "\t", row.names=1)
rownames(ethiopia.annot) <- gsub("-", ".", rownames(ethiopia.annot))
sudan.annot <- read.table("DNA/SihamAnalysis/Sample_SNP_location_fixed.txt", sep = "\t", row.names=1, header=TRUE)
rownames(sudan.annot) <- gsub(" ", "", rownames(sudan.annot))
rownames(sudan.annot) <- gsub("(", "", rownames(sudan.annot), fixed = TRUE)
rownames(sudan.annot) <- gsub(")", "", rownames(sudan.annot), fixed = TRUE)
#Swiss
uganda.annot <- read.table("Uganda/breeds.txt", sep = "\t", row.names=1)

annot.samples <- read.table(file = "Africa.Swiss.Adaptmap.samples.annotation.txt", sep = "\t", row.names=1, header=TRUE)
annot.samples <- cbind(annot.samples, Breed = NA, Continent = NA, Country = NA)

for(x in 1:nrow(annot.samples)){
  name = rownames(annot.samples)[x]
  if(name %in% rownames(algeria.annot)){
    idx <- which(rownames(algeria.annot) == name)
    annot.samples[x, "Breed"] <- as.character(algeria.annot[idx,"Race"])
    annot.samples[x, "Continent"] <- "Africa"
    annot.samples[x, "Country"] <- "Algeria"
  }
  if(name %in% rownames(ethiopia.annot)){
    idx <- which(rownames(ethiopia.annot) == name)
    annot.samples[x, "Breed"] <- as.character(ethiopia.annot[idx,1])
    annot.samples[x, "Continent"] <- "Africa"
    if(as.character(ethiopia.annot[idx,1]) == "Central Highland"){
      annot.samples[x, "Country"] <- "Cameroon"
    }else{
      annot.samples[x, "Country"] <- "Ethiopia"
    }
  }
  if(name %in% rownames(sudan.annot)){
    idx <- which(rownames(sudan.annot) == name)
    annot.samples[x, "Breed"] <- as.character(sudan.annot[idx,"Breed"])
    annot.samples[x, "Continent"] <- "Africa"
    annot.samples[x, "Country"] <- "Sudan"
  }
  if(annot.samples[x, 1] == "uganda"){
    breedID <- substr(name,1,3)
    idx <- which(rownames(uganda.annot) == breedID)
    annot.samples[x, "Breed"] <- as.character(uganda.annot[idx,1])
    annot.samples[x, "Continent"] <- "Africa"
    annot.samples[x, "Country"] <- "Uganda"
  }
  if(annot.samples[x, 1] == "adaptmap"){
    breedID <- unlist(strsplit(rownames(annot.samples)[x], "0"))[1]
    bIDs <- unlist(strsplit(breedID, "_"))
    
    idx <- which(rownames(adaptmap.annot) == paste0(bIDs[2], "_", bIDs[1]))
    if(length(idx) == 1){
      annot.samples[x, "Breed"] <- as.character(adaptmap.annot[idx,1])
      annot.samples[x, "Continent"] <- as.character(adaptmap.annot[idx,2])
      annot.samples[x, "Country"] <- as.character(adaptmap.annot[idx,3])
    }
  }
  if(annot.samples[x, 1] == "swiss"){
      annot.samples[x, "Continent"] <- "Europe"
      annot.samples[x, "Country"] <- "Switzerland"
  }
  if(annot.samples[x, 1] == "botswana"){
      annot.samples[x, "Continent"] <- "Africa"
      annot.samples[x, "Country"] <- "Botswana"
  }  
}

sudanIDX <- which(annot.samples[,1] == "sudan")
annot.samples[sudanIDX, 2] <- gsub("Nu", "Nubian", annot.samples[sudanIDX, 2])
annot.samples[sudanIDX, 2] <- gsub("Ni", "Nilotic", annot.samples[sudanIDX, 2])
annot.samples[sudanIDX, 2] <- gsub("Tagg", "Taggar", annot.samples[sudanIDX, 2])
annot.samples[sudanIDX, 2] <- gsub("Dese", "Desert goat", annot.samples[sudanIDX, 2])

annot.samples[sample(1:nrow(annot.samples), 100),]

write.table(annot.samples, file = "Africa.Swiss.Adaptmap.breed.annotation.txt", sep = "\t", quote=FALSE)
