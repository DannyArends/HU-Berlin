# pedigree.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014
#
# Analysis of genetic marker to discover wrong/imcompatible genotypes based on children and parents

setwd("D:/Stefan_Mouse_F3")
pedigreedata <- read.table("20140625_pedigree_fuerSebastian3.txt", sep="\t", header=TRUE)

OKorNotFinal <- NULL
for(r in 1:nrow(pedigreedata)){
  fatherRow <- which(pedigreedata[,"ID"] == pedigreedata[r, "Vater"])
  motherRow <- which(pedigreedata[,"ID"] == pedigreedata[r, "Mutter"])
  if(length(fatherRow) == 0) OKorNot <- rep(NA, 3)
  if(length(motherRow) == 0) OKorNot <- rep(NA, 3)
  if(length(fatherRow) == 1 && length(motherRow) == 1){
    cnt <- 1
    OKorNot <- rep(NA, 3)
    for(marker in c("ccna2", "bbs7", "m313")){
      ourMarker <- pedigreedata[r, marker]
      patMarker <- pedigreedata[fatherRow, marker]
      matMarker <- pedigreedata[motherRow, marker]
      if(is.na(ourMarker) || is.na(patMarker) || is.na(matMarker)){
        OKorNot[cnt] <- NA
      }else if(ourMarker == 0){
        if(patMarker == 0 && matMarker == 0) OKorNot[cnt] <- "000-0"
        if(patMarker == 0 && matMarker == 1) OKorNot[cnt] <- "001-0"
        if(patMarker == 0 && matMarker == 2) OKorNot[cnt] <- "002-1"        # 1 = maternal error
        if(patMarker == 1 && matMarker == 0) OKorNot[cnt] <- "010-0"
        if(patMarker == 2 && matMarker == 0) OKorNot[cnt] <- "020-2"        # 2 = paternal error
        if(patMarker == 1 && matMarker == 1) OKorNot[cnt] <- "011-0"
        if(patMarker == 1 && matMarker == 2) OKorNot[cnt] <- "012-3"
        if(patMarker == 2 && matMarker == 1) OKorNot[cnt] <- "021-3"        # 3 = parent error        
        if(patMarker == 2 && matMarker == 2) OKorNot[cnt] <- "022-3"
      }else if(ourMarker == 1){
        if(patMarker == 0 && matMarker == 0) OKorNot[cnt] <- "100-3"
        if(patMarker == 0 && matMarker == 1) OKorNot[cnt] <- "101-0"
        if(patMarker == 0 && matMarker == 2) OKorNot[cnt] <- "102-0"
        if(patMarker == 1 && matMarker == 0) OKorNot[cnt] <- "110-0"
        if(patMarker == 2 && matMarker == 0) OKorNot[cnt] <- "120-0"
        if(patMarker == 1 && matMarker == 1) OKorNot[cnt] <- "111-0"
        if(patMarker == 1 && matMarker == 2) OKorNot[cnt] <- "112-0"
        if(patMarker == 2 && matMarker == 1) OKorNot[cnt] <- "121-0"     
        if(patMarker == 2 && matMarker == 2) OKorNot[cnt] <- "122-3"
      }else if(ourMarker == 2){
        if(patMarker == 0 && matMarker == 0) OKorNot[cnt] <- "200-3"
        if(patMarker == 0 && matMarker == 1) OKorNot[cnt] <- "201-3"
        if(patMarker == 0 && matMarker == 2) OKorNot[cnt] <- "202-2"
        if(patMarker == 1 && matMarker == 0) OKorNot[cnt] <- "210-3"
        if(patMarker == 2 && matMarker == 0) OKorNot[cnt] <- "220-1"
        if(patMarker == 1 && matMarker == 1) OKorNot[cnt] <- "211-0"
        if(patMarker == 1 && matMarker == 2) OKorNot[cnt] <- "212-0"
        if(patMarker == 2 && matMarker == 1) OKorNot[cnt] <- "221-0"
        if(patMarker == 2 && matMarker == 2) OKorNot[cnt] <- "222-0"
      }
      cnt <- cnt + 1
    }
  }
  OKorNotFinal <- rbind(OKorNotFinal, OKorNot)
  cat(pedigreedata[r,"ID"], ":", OKorNot,"\n")    
}

rownames(OKorNotFinal) <- rownames(pedigreedata)
oldnames <- colnames(pedigreedata)
pedigreedata <- cbind(pedigreedata, OKorNotFinal)
colnames(pedigreedata) <- c(oldnames, "Rccna2","Rbbs7","Rm313")

Rccna2 <- as.numeric(unlist(lapply(strsplit(as.character(pedigreedata[,"Rccna2"]),"-"),"[",2)))
pedigreedata[which(Rccna2 > 0), c("ID","Wurf_ID","Rccna2")]

Rbbs7 <- as.numeric(unlist(lapply(strsplit(as.character(pedigreedata[,"Rbbs7"]),"-"),"[",2)))
pedigreedata[which(Rbbs7 > 0), c("ID","Wurf_ID","Rbbs7")]

Rm313 <- as.numeric(unlist(lapply(strsplit(as.character(pedigreedata[,"Rm313"]),"-"),"[",2)))
pedigreedata[which(Rm313 > 0), c("ID","Wurf_ID","Rm313")]