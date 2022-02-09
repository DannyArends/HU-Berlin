#
# Analysis Siham
#

setwd("D:/Edrive/Goat/KASP")
french <- read.table("French_All_Chr6.txt",sep="\t", header=TRUE, row.names=2)
vargoat <- read.table("VarGoat_All_Chr6.txt",sep="\t", header=TRUE, row.names=2)

# All french samples are in the vargoat, so we only need vargoat
all(colnames(french) %in% colnames(vargoat))

CSN1S1 <- read.table("CSN1S1_SNPtoPV_All.txt",sep="\t", header= TRUE, check.names=FALSE)
CSN2 <- read.table("CSN2_SNPtoPV_All.txt",sep="\t", header= TRUE, check.names=FALSE)
CSN1S2 <- read.table("CSN1S2_SNPtoPV_All.txt",sep="\t", header= TRUE, check.names=FALSE)
CSN3 <- read.table("CSN3_SNPtoPV_All.txt",sep="\t", header= TRUE, check.names=FALSE)

snps <- c(colnames(CSN1S1)[1:4], colnames(CSN2)[1], colnames(CSN1S2)[1:2], colnames(CSN3)[1:3])

for(x in 9:ncol(vargoat)){
  vargoat[snps,x] <- unlist(lapply(strsplit(vargoat[snps,x], ":"), "[",1))
}
vargoat <- vargoat[snps,]
vargoat[vargoat == "./."] <- NA

for(snp in rownames(vargoat)){
  hR <- paste0(vargoat[snp, "REF"], vargoat[snp, "REF"])
  hH <- "H"
  hA <- paste0(vargoat[snp, "ALT"], vargoat[snp, "ALT"])
  vargoat[snp, which(vargoat[snp, ] == "0/0")] <- hR
  vargoat[snp, which(vargoat[snp, ] == "0/1")] <- hH
  vargoat[snp, which(vargoat[snp, ] == "1/0")] <- hH
  vargoat[snp, which(vargoat[snp, ] == "1/1")] <- hA
}

vargoat <- vargoat[, -c(5:8)]

res <- matrix(NA, length(5:ncol(vargoat)), 4, dimnames=list(colnames(vargoat[,5:ncol(vargoat)]), c("CSN1S1", "CSN2", "CSN1S2", "CSN3")))
for(x in 5:ncol(vargoat)){
  gts <- vargoat[colnames(CSN1S1)[1:4],x]
  if(!any(is.na(gts))) {
    for(y in 1:nrow(CSN1S1)){
      if(all(CSN1S1[y,colnames(CSN1S1)[1:4]] == gts)){ 
        res[colnames(vargoat)[x], "CSN1S1"] <- CSN1S1[y, "Proteinvariant"]
      }
    }
  }
}

for(x in 5:ncol(vargoat)){
  gts <- vargoat[colnames(CSN2)[1],x]
  if(!any(is.na(gts))) {
    for(y in 1:nrow(CSN2)){
      if(all(CSN2[y,colnames(CSN2)[1]] == gts)){ 
        res[colnames(vargoat)[x], "CSN2"] <- CSN2[y, "Proteinvariant"]
      }
    }
  }
}

for(x in 5:ncol(vargoat)){
  gts <- vargoat[colnames(CSN1S2)[1:2],x]
  if(!any(is.na(gts))) {
    for(y in 1:nrow(CSN1S2)){
      if(all(CSN1S2[y,colnames(CSN1S2)[1:2]] == gts)){ 
        res[colnames(vargoat)[x], "CSN1S2"] <- CSN1S2[y, "Proteinvariant"]
      }
    }
  }
}

for(x in 5:ncol(vargoat)){
  gts <- vargoat[colnames(CSN3)[1:3],x]
  if(!any(is.na(gts))) {
    for(y in 1:nrow(CSN3)){
      if(all(CSN3[y,colnames(CSN3)[1:3]] == gts)){ 
        res[colnames(vargoat)[x], "CSN3"] <- CSN3[y, "Proteinvariant"]
      }
    }
  }
}

breeddata <- read.table("breeds.txt", sep = "\t", header=TRUE)


siham <- read.table("KASP_Lab_Siham.txt",sep = "\t", header=TRUE, na.strings=c("", "NA", "-"), row.names=1)
siham <- siham[!is.na(siham[,"AlphaS1casein"]),]

proteinsVars <- unique(c(siham[, "AlphaS1casein"], res[, "CSN1S1"]))


summarySiham <- matrix(NA, length(unique(siham[, "Breed"])), length(proteinsVars), dimnames=list(unique(siham[, "Breed"]), proteinsVars))

for(breed in unique(siham[, "Breed"])){
  bt <- table(siham[siham[, "Breed"] == breed,"AlphaS1casein"])
  summarySiham[breed, names(bt)] <- bt
}

breeddata <- cbind(mL = paste0(breeddata[,"Breed_Country.code"], "CH"),breeddata)

ss <- strsplit(rownames(res), ".", fixed=TRUE)
codez <- paste0(substr(unlist(lapply(ss,"[",2)),0,3), "_", unlist(lapply(ss,"[",1)))
idx <- which(codez %in% breeddata[, "mL"])
res <- res[idx, ]
rownames(res) <- codez[idx]

summaryvarGoat <- matrix(NA, length(unique(codez[idx])),  length(proteinsVars), dimnames=list(unique(codez[idx]), proteinsVars))
for(breed in unique(codez[idx])){
  bt <- table(res[which(rownames(res) == breed),"CSN1S1"])
  summaryvarGoat[breed, names(bt)] <- bt
}

mm <- rbind(summarySiham, summaryvarGoat)
mm <- mm[which(apply(mm,1,sum,na.rm=TRUE) >= 5),]
ma <- mm[,-which(apply(mm, 2, function(x){all(is.na(x))}))]

mP <- t(apply(ma,1, function(x){return(x/sum(x,na.rm=TRUE))}))
mP[is.na(mP)] <- 0
plot(hclust(dist(mP)))


alleles <- c("A","B1", "B2", "B3", "C", "C1", "J")

mAll <- matrix(0, nrow(ma), length(alleles), dimnames= list(rownames(ma), alleles))

for(x in 1:ncol(ma)){
  ss <- strsplit(colnames(ma)[x], "_")
  a1 <- ss[[1]][1]
  a2 <- ss[[1]][2]
  for(y in 1:nrow(ma)){
    if(!is.na(ma[y,x])) mAll[y, a1] <- mAll[y, a1] + ma[y,x]
    if(!is.na(ma[y,x])) mAll[y, a2] <- mAll[y, a2] + ma[y,x]
  }
}

mP <- t(apply(mAll,1, function(x){return(x/sum(x,na.rm=TRUE))}))

for(x in 7:nrow(mP)){
  idxBD <- which(breeddata[,1] == rownames(mP)[x])
  rownames(mP)[x] <- paste0(breeddata[idxBD, "Breed.name"], " (",breeddata[idxBD, "Country"], ")")
}

#plot(hclust(dist(mP)))
plot(hclust(dist(mP)), main = "CSN1S1 (breeds >= 5 animals)", xlab="Breed")

