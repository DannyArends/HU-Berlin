# Data Pre-processing
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written March, 2009
#

setwd("E:/Mouse/ClassicalPhenotypes/FV3")

crossFF <- read.table("cross_F2_FF_input.csv", header=TRUE, sep=",", colClasses="character")
crossFF[which(crossFF[,"Sex"]==0),"Sex"] <- 'f'
crossFF[which(crossFF[,"Sex"]==1),"Sex"] <- 'm'
crossFF[3:nrow(crossFF), grep("Mit", colnames(crossFF))] <- apply(crossFF[3:nrow(crossFF), grep("Mit", colnames(crossFF))],2,function(x){
  r <- rep(NA, length(x))
  r[which(x==1)] <- "A"
  r[which(x==2)] <- "H"
  r[which(x==3)] <- "B"
  r
})
write.table(crossFF,"cross_F2_FF.csv", quote=FALSE, row.names=FALSE, sep=",")

crossNF <- read.table(file="cross_F2_NF_input.csv", header=TRUE, sep=",", colClasses="character")
crossNF[3:nrow(crossNF), 63:ncol(crossNF)] <- apply(crossNF[3:nrow(crossNF), 63:ncol(crossNF)],2,function(x){
  r <- rep(NA, length(x))
  r[which(x=="1")] <- "A"
  r[which(x=="3")] <- "H"
  r[which(x=="2")] <- "B"
  r
})
write.table(crossNF,"cross_F2_NF.csv", row.names=FALSE, quote=FALSE, sep=",")

lociFF <- t(crossFF[1:2,27:ncol(crossFF)]) ; lociNF <- t(crossNF[1:2,63:ncol(crossNF)])                             # Get both the marker sets
loci <- rbind(lociFF,lociNF)
loci <- loci[unique(rownames(loci)),]                                                                               # Get the unique loci
colnames(loci) <- c("Chr", "Loc")                                                                                   # Add column identifiers

loci[loci[,"Chr"] == "X","Chr"] <- 20                                                                               # Rename X
loci <- loci[order(as.numeric(loci[,"Chr"]), as.numeric(loci[,"Loc"])),]                                            # Sort on Chr, Loc
loci[loci[,"Chr"] == "20","Chr"] <- "X"

IndFF <- crossFF[,"ID"][3:nrow(crossFF)]                                                                            # Individuals in FF
IndNF <- crossNF[,"ID"][3:nrow(crossNF)]                                                                            # Individuals in NF

genotypeMatrix <- matrix(NA,nrow(loci),length(IndFF) + length(IndNF))
genotypes <- cbind(loci, genotypeMatrix)
colnames(genotypes) <- c("Chr","Loc", IndFF, IndNF)                                                                 # We now have the empty combined Genotypes Matrix

for(x in IndFF){                                                                                                    # Fill in the FF animals
  markerdata  <- crossFF[crossFF[,"ID"] == x,27:ncol(crossFF)]
  markernames <- names(markerdata)
  genotypes[markernames,x] <- as.character(markerdata)
}

for(x in IndNF){                                                                                                    # Fill in the NF animals
  markerdata  <- crossNF[crossNF[,"ID"] == x,63:ncol(crossNF)]
  markernames <- names(markerdata)
  genotypes[markernames,x] <- as.character(markerdata)
}

genotypes <- t(genotypes)

crossFFm <- crossFF[,which(colnames(crossFF) %in% colnames(crossNF)[2:63])]
crossNFm <- crossNF[3:nrow(crossNF),match(colnames(crossFFm),colnames(crossNF))]

crossALL <- cbind(ID=c("","", rownames(genotypes)[-c(1:2)]), rbind(crossFFm,crossNFm), Futter = c("","", rep("0", length(IndFF)), rep("1", length(IndNF))), genotypes)
write.table(crossALL, "cross_F2.csv", row.names=FALSE, quote=FALSE, sep=",")

