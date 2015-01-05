# regionPlot.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of candidate gene data from Sebastiaan

setwd("E:/Mouse/DNA/DiversityArray/")
snpOUT  <- read.table("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
regions <- read.table("Analysis/Diabetes/BFMI861-S1andBFMI860-12vsALL_Regions.txt", sep="\t", header=TRUE)


asSNPs <- matrix(NA,nrow(snpOUT),length(9:14))
for(x in 1:nrow(snpOUT)){
  for(i in 9:14){
    if(is.na(snpOUT[x,i])){
      asSNPs[x,(i-8)] <- NA
    }else{
      if(snpOUT[x,i] == 0) asSNPs[x,(i-8)] <- paste0(snpOUT[x,"JAX_A"],snpOUT[x,"JAX_A"])
      if(snpOUT[x,i] == 1) asSNPs[x,(i-8)] <- paste0(unlist(sort(as.character(c(snpOUT[x,"JAX_A"],snpOUT[x,"JAX_B"])))), collapse="")
      if(snpOUT[x,i] == 2) asSNPs[x,(i-8)] <- paste0(snpOUT[x,"JAX_B"],snpOUT[x,"JAX_B"])
    }
  }
}

lvls <- levels(as.factor(asSNPs))    # "AA" "AC" "AG" "AT" "CC" "CG" "CT" "GG" "GT" "TT"
iupaccol <- c("green", "gray", "gray", "gray", "blue", "gray","gray", "black", "gray", "red")
numgeno <- NULL
for(x in 1:nrow(asSNPs)){
  numgeno <- rbind(numgeno, match(asSNPs[x,],lvls))
}

rownames(numgeno) <- snpOUT[,"JAX_ID"]
colnames(numgeno) <- colnames(snpOUT[9:14])

for(x in 1:nrow(regions)){
  png(paste0("Analysis/Figures/region",x,".png"), width=2000, height = 800)
  inRegion <- which(snpOUT[,"Chr"] == regions[x,"Chr"] & as.numeric(snpOUT[,"Location"]) > (regions[x,"Start"]-4000000) & as.numeric(snpOUT[,"Location"]) < (regions[x,"End"] + 4000000))
  ngeno <- numgeno[inRegion, ]
  op <- par(mar = c(5,  7, 4, 2) + 0.1)
  image(1:nrow(ngeno), 1:ncol(ngeno), ngeno, xlab="JAX", ylab="", yaxt="n", xaxt="n", col=iupaccol, main=paste0("Region: chr",regions[x,"Chr"],":", regions[x,"Start"],"-",regions[x,"End"]))
  axis(2, at=1:ncol(ngeno), colnames(ngeno),las=2, cex.axis=1.2)
  axis(1, at=1:nrow(ngeno), gsub("JAX00","", rownames(ngeno)), las=2, cex.axis=0.5)
  box()
  dev.off()
}

