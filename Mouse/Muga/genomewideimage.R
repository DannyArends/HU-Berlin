# genomewideimage.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014


setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                              # Normal A, H, B genotypes

BFMI <- genotypes[,"BFMI860-12 (V2)"]
names(BFMI) <- rownames(genotypes)

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                                                                    # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("Vater", "W.dat", "W.Label", "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
            
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)), phenos]                  # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F2 <- F2[-which(F2=="6661459")]                                                                               # This individual has no genotype data

onegenotype <- which(lapply(apply(genotypes[,F2], 1, table), length) == 1)                                    # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- unique(genotypes[-onegenotype, F2])                                                            # Only take the F2 individuals, and the unique markers
BFMI <- BFMI[rownames(genotypes)]
unabletodo <- which(is.na(BFMI))
unabletodo <- c(unabletodo, which(BFMI == "H"))

BFMI <- BFMI[-unabletodo]
genotypes <- genotypes[-unabletodo, ]
enough <- apply(genotypes,1,function(x){
  tbl <- na.omit(table(as.character(x)))
 return((!any(tbl < 10)) && length(tbl) >= 2)
})
genotypes <- genotypes[which(enough),]
BFMI <- BFMI[which(enough)]
cat("Left with", nrow(genotypes),",BFMI:", length(BFMI), "markers\n")                                                  # == Left with 11677 markers


nind <- length(genotypes) + 1
numericgenotypes <- t(apply(cbind(BFMI, genotypes), 1, function(x){ return(as.numeric(x[2:nind] == x[1])) }))
colnames(numericgenotypes) <- colnames(genotypes)

chr3marker <- "UNC5048297"
plot(phenotypes[F2,"mri70d_fat"] ~ as.factor(unlist(genotypes[chr3marker,]))) # "A"
heavy <- names(which(unlist(genotypes[chr3marker,]) == "A"))
other <- names(which(unlist(genotypes[chr3marker,]) != "A"))

numericgenotypes[which(genotypes == "H")] <- 0.5

mlabels <-  apply(map[rownames(genotypes),c(1,2)], 1, function(x){round(as.numeric(x[2])/1000000,3) })

heavyS <- sort(phenotypes[colnames(genotypes[,heavy]),"mri70d_fat"],index.return=TRUE)
otherS <- sort(phenotypes[colnames(genotypes[,other]),"mri70d_fat"],index.return=TRUE)

rmap <- map[rownames(numericgenotypes),]

for(chr in unique(rmap[rownames(numericgenotypes),"Chr"])){
  onChr <- rownames(rmap[which(rmap[,"Chr"] == chr),])
  ngeno <- numericgenotypes[onChr,]
  if(nrow(ngeno) >= 2){
    png(paste0("Sorted_Chr",chr,".png"), width=max(nrow(ngeno)*20,800), height = nind*20)
      op <- par(mar = c(5, 8, 4, 2) + 0.1)
      image(x=1:nrow(ngeno), y=1:nind, ngeno[, c(other[otherS$ix], heavy[heavyS$ix])], yaxt = "n", xaxt = "n", col=c("yellow", "gray","red"), ylab="", xlab="")
      abline(v = 1:nrow(ngeno) + 0.5, col="white", lwd = 0.5)
      abline(h = 1:nind + 1, col="white", lwd = 0.5)
      abline(h = length(other)+1, col="black", lwd = 5)
      axis(1, at = 1:nrow(ngeno), as.character(mlabels[onChr]), las=2, cex.axis=0.9)
      axis(2, at = (1:(nind-1) + 0.5), phenotypes[c(other[otherS$ix], heavy[heavyS$ix]),"Vater"], las=2, cex.axis=0.9)
      box()
    dev.off()
  }
}

