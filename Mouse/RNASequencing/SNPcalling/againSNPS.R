
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")

getDP <- function(x){
  v1 <- strsplit(x, ";")
  as.numeric(unlist(strsplit(gsub("DP=","",unlist(v1)[grep("DP=", unlist(v1))]),",")))
}

createNames <- function(x){ paste0(x[,1],":", x[,2],"_", x[,5]) }
getGenotypes <- function(x){ unlist(lapply(strsplit(x, ":"),"[",1)) }
getProbabilities <- function(x){ unlist(lapply(strsplit(x, ":"),"[",2)) }
getMaxProb <- function(x){ 
  p <- getProbabilities(x)
  unlist(lapply(strsplit(p, ","),function(i){max(as.numeric(i))}))
}

##### Analysis of the SNPs

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
vcfdata <- read.table("population.vcf", colClasses="character")
colnames(vcfdata) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","5068","5069","4868","5067","5070","5071","5072","5073","5074","5075")
rownames(vcfdata) <- createNames(vcfdata)
parents <- c("5068","5069","4868","5067")
matB6N  <- c("5070","5071","5072")
matBFMI <- c("5073","5074","5075")
samples <- c("5068","5069","4868","5067","5070","5071","5072","5073","5074","5075")

DPs <- unlist(lapply(vcfdata[,"INFO"], getDP))
vcfdata <- vcfdata[DPs > 100,]

onAutosomes <- length(which(!(vcfdata[,"CHROM"] == "X" | vcfdata[,"CHROM"] == "Y" | vcfdata[,"CHROM"] == "MT")))
cat("Called",onAutosomes,"/",dim(vcfdata)[1],"variants (autosomes / all)\n")

genos <- matrix(NA, nrow(vcfdata), length(samples))
probs <- matrix(NA, nrow(vcfdata), length(samples))
colnames(genos) <- samples ; rownames(genos) <- createNames(vcfdata)
colnames(probs) <- samples ; rownames(probs) <- createNames(vcfdata)
for(s in samples){ 
  genos[,s] <- getGenotypes(vcfdata[,s])
  probs[,s] <- getMaxProb(vcfdata[,s])
}


confidenceThreshold <- 50

probInParents <- which(apply(probs[,parents], 1,function(x){ sum(x > confidenceThreshold) == 4 }))
genos <- genos[probInParents, ] ; probs <- probs[probInParents, ]

chroms <- unlist(lapply(strsplit(rownames(genos),":"),"[",1))
onAutosomes <- length(which(!(chroms == "X" | chroms == "Y" | chroms == "MT")))

cat("Detected",onAutosomes,"/",dim(genos)[1],"variants\n")

difGenoInParents <- which(apply(genos[,parents], 1,function(x){ 
  if(any(x == "0/1")) return(FALSE)
  if(x[1] == x[2] && x[3] == x[4] && x[1] != x[3]) return(TRUE)
  return(FALSE)
}))

genos <- genos[difGenoInParents, ] ; probs <- probs[difGenoInParents, ]
chroms <- unlist(lapply(strsplit(rownames(genos),":"),"[",1))
onAutosomes <- length(which(!(chroms == "X" | chroms == "Y" | chroms == "MT")))
cat("Detected",onAutosomes,"/",dim(genos)[1],"variants usable for ASE\n")
