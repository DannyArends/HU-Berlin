setwd("D:/Edrive/Mouse/MDC/Jan2020")

geno <- read.csv("7984_geno.txt", sep = "\t")
pheno <- read.csv("7984_pheno.txt", sep = "\t", na.strings = c("", "NA", "-", "X"))

colnames(geno) <- c("MouseID", "GT", "Uncertain")

mother <- NA
father <- NA
wdate <- NA
WG <- NA
GEN <- NA
for(x in 1:nrow(pheno)){
  if(!is.na(pheno[x, "Vater"])){
    mother <- pheno[x, "Mutter"]
    father <- pheno[x, "Vater"]
    wdate <- pheno[x, "W.dat"]
    WG <- pheno[x, "WG"]
    GEN <- pheno[x, "Gen."]
  }else{
    pheno[x, "Mutter"] <- mother
    pheno[x, "Vater"] <- father
    pheno[x, "W.dat"] <- wdate
    pheno[x, "WG"] <- WG
    pheno[x, "Gen."] <- GEN
  }
}

sex <- as.character(pheno[,"m"])
sex[which(is.na(sex))] <- "f"
pheno <- cbind(pheno, sex = sex)

pheno <- pheno[, c(1:5, 10:18,25)]

pheno <- cbind(pheno, GenGood = NA)
lookupTable <- c()
for(x in 1:nrow(pheno)){
  mother <- gsub("MDC-", "", as.character(pheno[x, "Mutter"]))
  father <- gsub("MDC-", "", as.character(pheno[x, "Vater"]))
  child <- gsub("MDC-", "", as.character(pheno[x, "ID.Nr"]))
  
  mGen <- NA
  fGen <- NA
  cGen <- NA

  if(substr(mother, 0, 3) == "101") mGen <- 0
  if(substr(mother, 0, 3) == "TCF") mGen <- 0
  if(substr(mother, 0, 3) == "860") mGen <- 0
  
  if(substr(father, 0, 3) == "101") fGen <- 0
  if(substr(father, 0, 3) == "TCF") fGen <- 0
  if(substr(father, 0, 3) == "860") fGen <- 0
  
  if(is.na(mGen)){
    #Figure out what generation the mother is
    idx <- which(lookupTable[,1] == mother)
    if(length(idx == 1)){
      mGen <- as.numeric(lookupTable[idx, 2])
    }else{
      stop(paste0("mother ",mother," not in lookup table"))
    }
  }
  if(is.na(fGen)){
    #Figure out what generation the father is
    idx <- which(lookupTable[,1] == father)
    if(length(idx == 1)){
      fGen <- as.numeric(lookupTable[idx, 2])
    }else{
      stop(paste0("father ",father," not in lookup table"))
    }
  }
  if(!is.na(mGen) && !is.na(fGen) && mGen == 0 && fGen == 0){
    cGen <- 1
    lookupTable <- rbind(lookupTable, c(child, cGen))
  }
  if(!is.na(mGen) && !is.na(fGen) && max(mGen, fGen) > 0){
    cGen <- max(mGen, fGen) + 1
    lookupTable <- rbind(lookupTable, c(child, cGen))
  }
  pheno[x, "GenGood"] <- cGen
}

# set the name of the animal as the rowname !!!
rownames(pheno) <- gsub("MDC-", "", pheno[, "ID.Nr"])
rownames(geno) <- paste0("7984-", geno[, "MouseID"])

mdata <- cbind(pheno, geno[rownames(pheno),])

mdata <- mdata[-which(is.na(mdata[, "GT"])),]
mdata <- mdata[-which(!is.na(mdata[, "Uncertain"])),]
mdata <- mdata[, -ncol(mdata)]

mdata <- mdata[-which(apply(apply(mdata, 1, is.na),2,sum) > 0),]

idxs <- which(mdata[, "GenGood"] == 10)
boxplot(mdata[idxs, "X70"] ~ mdata[idxs, "GT"])