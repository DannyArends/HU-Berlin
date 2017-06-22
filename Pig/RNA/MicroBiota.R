setwd("D:/Edrive/Pig/RNA/Zink_mRNASequencing")
rpkm <- read.table("RPKMnorm.txt", sep="\t", check.names=FALSE)
colnames(rpkm) <- unlist(lapply(strsplit(colnames(rpkm),".",fixed=TRUE), "[",1))

# Groups
high <- c("92", "113", "116")

rpkm[rpkm == 0] <- NA
rpkm <- log2(rpkm)
rpkm[is.na(rpkm)] <- 0

cat("Before:", nrow(rpkm), "genes\n")
rpkm <- rpkm[-which(apply(rpkm[,medium],1,function(x){all(x < 1)})),]
rpkm <- rpkm[-which(apply(rpkm[,high],1,function(x){all(x < 1)})),]

cat("After:", nrow(rpkm), "genes\n")

setwd("D:/Edrive/Pig/Microbiota")

microbiota <- read.table("totalreads.txt",sep="\t",header=TRUE)

microbiota <- microbiota[which(microbiota[,"Tier"] %in% c(medium, high)),]
op <- par(mfrow = c(1,3))
for(x in c(high)){
  cat("animal: ", x, "\n")
  inTier <- microbiota[which(microbiota[,"Tier"] == x),]
  totalAbundance <- sum(inTier[,"abundance"])
  inTier <- cbind(inTier, relAbundance = inTier[,"abundance"] / totalAbundance)
  inTier <- cbind(inTier, genusAbundance = NA)
  inTier <- cbind(inTier, genusRelAbundance = NA)
  genusAbundance <- NULL
  for(y in unique(inTier[, "genus"])){
    inGenus <- which(inTier[, "genus"] == y)
    inTier[inGenus, "genusAbundance"] <- sum(inTier[inGenus, "abundance"])
    inTier[inGenus, "genusRelAbundance"] <- sum(inTier[inGenus, "abundance"])/ totalAbundance
    genusAbundance <- c(genusAbundance, sum(inTier[inGenus, "abundance"])/ totalAbundance)
  }
  names(genusAbundance) <- unique(inTier[, "genus"])
  lowRepresentation <- which(genusAbundance < 0.01)
  tooLow <- sum(genusAbundance[lowRepresentation])
  genusAbundance <- genusAbundance[-lowRepresentation]
  genusAbundance <- c(genusAbundance, lowRepresentation = tooLow)
  pie(genusAbundance)
}

