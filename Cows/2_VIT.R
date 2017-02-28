
### Load in the files from the VIT

genotypes <- read.table("VIT/gZWSgenotypenDSNbullen_1612", colClasses="character",stringsAsFactor=FALSE)
gts <- apply(genotypes, 1, function(x){
  unlist(strsplit(x[2], ""))
})
colnames(gts) <- genotypes[,1]
vit_info <- read.table("VIT/SNPNamenV2mit_extra_Info_repaired.txt", sep=";", colClasses="character",stringsAsFactor=FALSE, header=TRUE)
vit_info <- vit_info[order(as.numeric(vit_info[,"SNP.Index"])), ]
rownames(gts) <- vit_info[,"SNP.Name"]
rownames(vit_info) <- vit_info[,"SNP.Name"]

for(x in 1:nrow(gts)){
  if(vit_info[rownames(gts)[x],"ALLELEA_TOP"] != ""){
    cod2 <- paste0(vit_info[rownames(gts)[x],"ALLELEA_TOP"], vit_info[rownames(gts)[x],"ALLELEA_TOP"])
    cod1 <- paste0(sort(c(vit_info[rownames(gts)[x],"ALLELEA_TOP"], vit_info[rownames(gts)[x],"ALLELEB_TOP"])), collapse="")
    cod0 <- paste0(vit_info[rownames(gts)[x],"ALLELEB_TOP"], vit_info[rownames(gts)[x],"ALLELEB_TOP"])
    gts[x, gts[x,] == "0"] <- rep(cod0, length(which(gts[x,] == "0")))
    gts[x, gts[x,] == "1"] <- rep(cod1, length(which(gts[x,] == "1")))
    gts[x, gts[x,] == "2"] <- rep(cod2, length(which(gts[x,] == "2")))
    gts[x, gts[x,] == "9"] <- rep(NA, length(which(gts[x,] == "9")))
  }else{
    gts[x,] <- rep(NA, ncol(gts))
  }
}
write.table(gts, file="GT_VIT.txt", sep="\t", quote=FALSE)

