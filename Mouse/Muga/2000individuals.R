setwd("E:/Mouse/DNA/MegaMuga")
chr12 <- read.table("chromsome12.txt", sep="\t", header=TRUE, check.names=FALSE,na.strings=c("","-"))
map <- t(chr12[1:2,])[-1,]
chr12 <- chr12[-c(1,2),]

individualsHaveOne <- as.character(chr12[which(apply(apply(chr12[,2:4],1,is.na),2,sum) != 3),1])
individualsHaveOne <- individualsHaveOne[-which(grepl("B", individualsHaveOne))]

mapO <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypesO   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                                    # Normal A, H, B genotypes

genotypesOG <- apply(genotypesO, 1, function(x){
  cat(x["BFMI860-12 (V2)"], x["B6N"],"\n")
  a <- rep("",length(x))
  a[as.character(x) == as.character(x["BFMI860-12 (V2)"])] <- "BFMI"
  a[as.character(x) != as.character(x["BFMI860-12 (V2)"])] <- "B6N"
  a[as.character(x) == "H"] <- "H"
  return(a)
})

chr12a <- apply(chr12[,2:4], 2, function(x){
  cat(x["474"], x["482"],"\n")
  a <- rep("",length(x))
  a[as.character(x) == as.character(x["474"])] <- "BFMI"
  a[as.character(x) != as.character(x["474"])] <- "B6N"
  a[as.character(x) == "H"] <- "H"
  return(a)
})

chr12a <- cbind(as.character(chr12[,1]), chr12a)

rownames(genotypesOG) <- colnames(genotypesO)

strt <-  82085109
stp <-  strt + 18891872

inR <- rownames(mapO[which(mapO[,1] == 12 & as.numeric(mapO[,2]) > strt & as.numeric(mapO[,2]) < stp),])

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
all2260 <- read.table("Data.txt",sep="\t",header=TRUE)

qtl <- read.table("Analysis/qtls_phased_RF1.txt",   sep="\t")
haveScore <- rownames(qtl[inR,][which(!is.na(qtl[inR, "marker"])),])

genotypesN <- chr12a[which(chr12a[,1] %in% individualsHaveOne),]
rownames(genotypesN) <- genotypesN[,1]
genotypesN <- cbind(map, t(genotypesN[,-1]))
colnames(genotypesN)[1:2] <- c("Chr","Mb_NCBI38")

genotypesOO <- genotypesOG[individualsHaveOne,haveScore]
genotypesOO <- cbind(mapO[colnames(genotypesOO),c(1,2)],t(genotypesOO))

#Bind together

combined <- rbind(genotypesOO,genotypesN)
combined <- combined[sort(as.numeric(combined[,2]),index.return=TRUE)$ix,]
combined <- cbind(LOD = qtl[rownames(combined),"marker"],combined)


write.table(combined,file="GenotypesChr12.txt",sep="\t")

phased <- genotypesPh[rownames(combined),]
phased <- phased[,-which(apply(phased, 2, function(x){return(sum(is.na(x)))}) > 70)]

phased <- cbind(LOD = qtl[rownames(phased),"marker"],phased)
write.table(phased,file="GenotypesChr12Phased.txt",sep="\t")