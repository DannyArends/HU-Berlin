# Combined all the different sources of horse data together for analysis of the horse SNP chip data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Nov, 2016
# first written Feb, 2015

toNumeric <- function(genotypes){
  numericG <- apply(genotypes, 1, function(x){
    geno <- table(unlist(strsplit(as.character(x),"")))
    a1 <- paste0(names(geno)[1],names(geno)[1])
    a2a <- paste0(c(names(geno)[1],names(geno)[2]),collapse="")
    a2b <- paste0(c(names(geno)[2],names(geno)[1]),collapse="")
    a3 <- paste0(names(geno)[2],names(geno)[2])
    ngeno <- rep(NA,length(x))
    ngeno[x == a1] <- 1
    ngeno[x == a2a] <- 2
    ngeno[x == a2b] <- 2
    ngeno[x == a3] <- 3
    return(ngeno)
  })
  rownames(numericG) <- colnames(genotypes)
  return(t(numericG))
}

toAB <- function(genotypes){
  numericG <- apply(genotypes, 1, function(x){
    geno <- table(unlist(strsplit(as.character(x),"")))
    a1 <- paste0(names(geno)[1],names(geno)[1])
    a2a <- paste0(c(names(geno)[1],names(geno)[2]),collapse="")
    a2b <- paste0(c(names(geno)[2],names(geno)[1]),collapse="")
    a3 <- paste0(names(geno)[2],names(geno)[2])
    ngeno <- rep(NA,length(x))
    ngeno[x == a1]  <- "AA"
    ngeno[x == a2a] <- "AB"
    ngeno[x == a2b] <- "AB"
    ngeno[x == a3]  <- "BB"
    return(ngeno)
  })
  rownames(numericG) <- colnames(genotypes)
  return(t(numericG))
}

setwd("E:/Horse/DNA/")

markerinfo  <- read.csv("60Karray/combinedmarkers.txt", sep="\t", na = "")
petersen    <- read.csv("Petersen2013/analysis/genotypes_snp.txt", sep="\t", na = c("NA", ""), row.names=NULL, colClasses="character")
arabian     <- read.csv("Equine60k/analysis/genotypes_snp.txt", sep="\t", na = c("NA", ""), row.names=NULL, colClasses="character")
kabadiner1  <- read.csv("Kabadiner/analysis/genotypes_snp_kabadiner-1.txt", sep="\t", na = c("NA", ""), row.names=NULL, colClasses="character")
kabadiner2  <- read.csv("Kabadiner/analysis/genotypes_snp_kabadiner-2.txt", sep="\t", na = c("NA", ""), row.names=NULL, colClasses="character")

all(markerinfo[,"Name"] == petersen[,1],na.rm=TRUE)   # Does the ordering match ?
all(petersen[,1] == arabian[,1],na.rm=TRUE)           # Does the ordering match ?
all(petersen[,1] == kabadiner1[,1],na.rm=TRUE)        # Does the ordering match ?
all(petersen[,1] == kabadiner2[,1],na.rm=TRUE)        # Does the ordering match ?

combined    <- cbind(petersen, arabian[,-1], kabadiner1[,-1], kabadiner2[,-1])

nolocation  <- which(is.na(markerinfo[,"MapInfo"]))
combined    <- combined[-nolocation, ]
rownames(combined) <- combined[,1]
combined <- combined[,-1]

markerinfo  <- markerinfo[-nolocation, ]
rownames(markerinfo) <- markerinfo[,"Name"]

cat("Starting with", nrow(combined), "probes\n")
nmissing <- apply(combined, 1, function(x){ return(sum(is.na(x))) })
enoughdata <- which((nmissing / ncol(combined)) <= 0.05)

combined <- combined[enoughdata, ]
markerinfo <- markerinfo[enoughdata, ]

cat(nrow(combined), "probes with 5% or less missing data\n")

noinfo  <- which(unlist(lapply(apply(combined, 1, table),length)) < 2)             # Non informative, since we only have 1 genotype

combined <- combined[-noinfo, ]
markerinfo <- markerinfo[-noinfo, ]

cat(nrow(combined), "probes are seggregating\n")

noduplicate   <- which(!duplicated(combined))                         # Duplicated markers
combined      <- combined[noduplicate, ]
markerinfo    <- markerinfo[noduplicate, ]

cat(nrow(combined), "probes are non-duplicates\n")

allelefreq <- apply(combined, 1 , function(x){
  tbl <- table(unlist(lapply(x, strsplit, "")))
  min(tbl / sum(tbl))
})

keep <- which(allelefreq >= 0.05)
combined     <- combined[keep, ]
markerinfo   <- markerinfo[keep, ]

cat(nrow(combined), "probes have an allele frequency above 0.05\n")

write.table(combined,   "combined/input/genotypes_snp.txt", sep="\t", quote=FALSE)
write.table(markerinfo, "combined/input/map.txt", sep="\t", quote=FALSE)
write.table(toNumeric(combined), file="combined/input/genotypes_num.txt", sep = "\t")
write.table(toAB(combined), file="combined/input/genotypes_AB.txt", sep = "\t")

h1   <- read.table("Kabadiner/input/kabadiner-1.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=2, nrow=15)
h2   <- read.table("Kabadiner/input/kabadiner-2.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=2, nrow=15)
a1   <- read.table("Equine60k/input/arabianhorses.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=1, nrow=28)
pe   <- read.table("Petersen2013/analysis/KSF1385391590.txt", sep="\t", header=FALSE)

a1 <- a1[,grep("GType", colnames(a1))]
colnames(a1) <- gsub(".GType", "", colnames(a1))
h1 <- h1[,grep("GType", colnames(h1))]
colnames(h1) <- gsub(".GType", "", colnames(h1))
h2 <- h2[,grep("GType", colnames(h2))]
colnames(h2) <- gsub(".GType", "", colnames(h2))
colnames(h2)[which(colnames(h2) == "P3611")] <- "P3611.1"

phenonames <- unique(c(rownames(a1),rownames(h1),rownames(h2)))

phenotypes <- matrix(NA, length(colnames(combined)), length(phenonames), dimnames=list(colnames(combined), phenonames))

for(x in 1:ncol(a1)){
  phenotypes[colnames(a1)[x], rownames(a1)] <- a1[,x]
}
for(x in 1:ncol(h1)){
  phenotypes[colnames(h1)[x], rownames(h1)] <- h1[,x]
}
for(x in 1:ncol(h2)){
  phenotypes[colnames(h2)[x], rownames(h2)] <- h2[,x]
}

for(x in 1:nrow(pe)){
  phenotypes[which(grepl(as.character(pe[x,1]), rownames(phenotypes))), "Strain"] <- as.character(pe[x,2])
}

i <- 778:849
phenotypes[i,"Strain"] <- phenotypes[i,"Breed"]

write.table(phenotypes, "combined/input/phenotypes.txt", sep="\t", quote=FALSE)
