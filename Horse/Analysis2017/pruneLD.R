library(genetics)
allelefreq <- apply(genotypes, 1 , function(x){
  tbl <- table(unlist(lapply(x, strsplit, "")))
  min(tbl / sum(tbl))
})

map.pruned <- cbind(map, "MAF" = allelefreq)
genotypes.pruned <- genotypes
markers <- rownames(genotypes.pruned)

x <- 1
while(x < length(markers)) {
  mName <- markers[x]
  mChr <- as.character(map[mName, "Chromosome"])
  mPos <- as.numeric(map[mName, "Position"])

  nearby <- which(as.character(map.pruned[,"Chromosome"]) == mChr & as.numeric(map.pruned[, "Position"]) > (mPos - 1000000) &  as.numeric(map.pruned[, "Position"]) < (mPos + 1000000))
  nearby <- rownames(map.pruned)[nearby]

  mGeno <- genotype(as.character(genotypes.pruned[mName,]), sep = "")
  locOfMarker <- which(nearby == mName)
  if(locOfMarker > 1) nearby <- nearby[-(1:locOfMarker-1)]

  LDs <- rep(NA, length(nearby))
  names(LDs) <- nearby
  for(y in nearby){
     LDs[y] <- round(LD(mGeno, genotype(as.character(genotypes.pruned[y,]), sep = ""))$"R^2",2)
  }
  inLD <- names(which(LDs > 0.5))
  if(length(inLD) > 1){
    MAFs <- map.pruned[inLD, ]
    bestSNP <- rownames(MAFs)[which.max(MAFs[, "MAF"])]
    inLD <- inLD[-which(inLD == bestSNP)]
    toPrune <- which(rownames(genotypes.pruned) %in% inLD)
    if(length(toPrune) > 0){
      genotypes.pruned <- genotypes.pruned[-toPrune, ]
      map.pruned <- map.pruned[-toPrune, ]
      markers <- rownames(genotypes.pruned)
    }
  }
  cat("Done", x, "/", nrow(genotypes), "/", nrow(genotypes.pruned), "==", nrow(map.pruned), "\n")
  x <- (x + 1)
}



### Removal based on LD
#library(genetics)
#genotypes[is.na(genotypes)] <- ""

#LDs <- matrix(NA, nrow(genotypes), nrow(genotypes), dimnames=list(rownames(genotypes), rownames(genotypes)))
#for(x in 1:nrow(genotypes)){
#  mChr <- as.character(map[x, "Chromosome"])
#  mPos <- as.numeric(map[x, "Position"])
#  nearby <- which(as.character(map[,"Chromosome"]) == mChr & as.numeric(map[, "Position"]) > (mPos - 1000000) &  as.numeric(map[, "Position"]) < (mPos + 1000000))
#  mGeno <- genotype(as.character(genotypes[x,]), sep = "")
#  nearby <- rownames(map)[nearby[-which(nearby <= x)]]
#  for(y in nearby){
#     LDs[x, y] <- round(LD(mGeno, genotype(as.character(genotypes[y,]), sep = ""))$"R^2",2)
#  }
#  if(x %% 100 == 0) cat("Done", x, "/", nrow(genotypes), "\n")
#}

#inLD <- c()
#for(x in 1:ncol(LDs)){
#  maxLD <- max(LDs[,x], na.rm=T)
#  if(maxLD >= 0.5) inLD <- c(inLD, x)
#}

#genotypes <- genotypes[-inLD,]
#map <- map[-inLD,]