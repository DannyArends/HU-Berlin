#
# Pre-process the horse data into a single format
#

opposite <- function(x){
  if(is.na(x)) return("")
  ret <- NULL
  for(e in x){
    if(e == "A") ret <- c(ret, "T")
    if(e == "C") ret <- c(ret, "G")
    if(e == "G") ret <- c(ret, "C")
    if(e == "T") ret <- c(ret, "A")
  }
  return(ret)
}

petersen_loc <- "Petersen2013/raw/"
setwd("E:/Horse/DNA/")

# Load raw data
fnames <- dir(petersen_loc)
files <- paste0(petersen_loc, fnames)

alldata <- vector("list", length(files))
names(alldata) <- fnames
for(x in 1:length(files)) {
  fp <- gzfile(files[x], "r")
  alldata[[fnames[x]]] <- read.csv(fp, header = FALSE)
  close(fp)
  cat("Loaded ", fnames[x],"\n")
}

# Load marker information (Illumina)
markerinfo <- read.csv("60Karray/GGP_Equine.csv", skip=7, header=TRUE, colClasses="character")
markerinfo[,"Name"] <- gsub("-", "_", as.character(markerinfo[,"Name"]))
cat("Loaded:", nrow(markerinfo), "markers\n")

# Load marker information from our own data
markerinfo2 <- read.csv("60Karray/equineSNP50.txt", header=TRUE, sep="\t", colClasses="character")
extra <- which(!(markerinfo2[,"Name"] %in% markerinfo[,"Name"]))
markerinfo2 <- markerinfo2[extra,]
extramarkers <- matrix(NA, length(extra), ncol(markerinfo))
colnames(extramarkers) <- colnames(markerinfo)
extramarkers[,"Name"]    <- markerinfo2[,"Name"]
extramarkers[,"Chr"]     <- markerinfo2[,"Chr"]
extramarkers[,"MapInfo"] <- markerinfo2[,"Position"]

cat("Adding", nrow(extramarkers), "markers with location information\n")
markerinfo <- rbind(markerinfo, extramarkers)

# Markers in Petersen's SNP data but not on the array
animals <- unique(unlist(lapply(alldata, function(x){ return(unique(x[,2])); })))
markersD <- gsub("-", "_", unique(unlist(lapply(alldata, function(x){ return(unique(x[,1])); }))))

extra <- which(!(markersD %in% markerinfo[,"Name"]))

extramarkers <- matrix(NA, length(extra), ncol(markerinfo))
colnames(extramarkers) <- colnames(markerinfo)
extramarkers[,"Name"]    <- markersD[extra]
cat("Adding", nrow(extramarkers), "markers without any information\n")
markerinfo <- rbind(markerinfo, extramarkers)

# Write out the combined marker information
write.table(markerinfo, "60Karray/combinedmarkers.txt", sep="\t", na = "")

markers <- markerinfo[,"Name"]

mdata <- matrix(NA, length(markers), length(animals), dimnames=list(markers, animals))
for(x in 1:length(alldata)){
  canimals <- unique(alldata[[x]][,2])
  for(anim in canimals){
    animsubset <- alldata[[x]][which(alldata[[x]][,2] == anim),]
    animsubset <- animsubset[which(animsubset[,7] >= 0.15),]
    mnames <- gsub("-", "_", as.character(animsubset[,1]))
    mdata[mnames, anim] <- apply(animsubset[,c(5,6)], 1, paste0, collapse="")
  }
  cat("Done",x," ", length(alldata),"\n")
}
write.table(mdata, "Petersen2013/analysis/genotypes_snp.txt", sep="\t", quote=FALSE, na = "")

# Our data (Arabian horses)

arabian <- read.table("Equine60k/input/arabianhorses.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=1)
genotypes             <- arabian[29:nrow(arabian), ]                  # Row 30 till the end contains genotype data
rownames(genotypes)   <- gsub("-", "_", arabian[29:nrow(arabian), 1])
genotypes <- genotypes[which(genotypes[,"GenTrain.Score"] >= 0.6), ]
genotypes <- genotypes[,c(grep("Top.Alleles", colnames(genotypes)),grep("TOP", colnames(genotypes)))]
colnames(genotypes) <- gsub(".TOP", "", colnames(genotypes))
colnames(genotypes) <- gsub(".Top.Alleles", "", colnames(genotypes))

arabGeno <- matrix(NA, nrow(mdata), ncol(genotypes), dimnames = list(rownames(mdata), colnames(genotypes)))

for(x in 1:nrow(arabGeno)){
  mname <- rownames(arabGeno)[x]
  genoIndex <- which(rownames(genotypes) == mname)
  cat(mname, genoIndex, "\n")
  if(length(genoIndex) == 1){
    genos <- unlist(genotypes[genoIndex, ])
    infoIndex <- which(markerinfo[,"Name"] == mname)
    if(length(infoIndex) == 1 && !is.na(markerinfo[infoIndex, "IlmnStrand"])){
      if(markerinfo[infoIndex, "IlmnStrand"] == "BOT" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        arabGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else if(markerinfo[infoIndex, "IlmnStrand"] == "TOP" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        arabGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else{
        arabGeno[mname, ] <- genos
      }
    }else{
      arabGeno[mname, ] <- genos
    }
  }
}
write.table(arabGeno, "Equine60k/analysis/genotypes_snp.txt", sep="\t", quote=FALSE, na = "")


refAlleles <- unlist(lapply(strsplit(markerinfo[,"SNP"],""),"[",2))
glabels <- apply(mdata,1,function(x){return(table(unlist(strsplit(x,""))))})

for(x in 1:length(glabels)){
  if(any(refAlleles[x] %in% names(glabels[[x]]))){
    cat("MATCH\n")
  }else{
    cat(markerinfo[x,"Name"], "No MATCH",refAlleles[x], "-> ",names(glabels[[x]]),"\n")
  }
}

# Our data (Kabadiner horses) 2

horsedata   <- read.table("Kabadiner/input/kabadiner-2.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=2)
genotypes   <- horsedata[15:nrow(horsedata), ]                                      # Row 14 till the end contains genotype data
rownames(genotypes)   <- gsub("-", "_", rownames(horsedata)[15:nrow(horsedata)])
genotypes <- genotypes[which(genotypes[,"Gentrain"] >= 0.6), ]
genotypes <- genotypes[,grep("TOP", colnames(genotypes))]
colnames(genotypes) <- gsub(".TOP", "", colnames(genotypes))

kabadinerGeno <- matrix(NA, length(markers), ncol(genotypes), dimnames = list(markers, colnames(genotypes)))

for(x in 1:nrow(kabadinerGeno)){
  mname <- rownames(kabadinerGeno)[x]
  genoIndex <- which(rownames(genotypes) == mname)
  cat(mname, genoIndex, "\n")
  if(length(genoIndex) == 1){
    genos <- unlist(genotypes[genoIndex, ])
    infoIndex <- which(markerinfo[,"Name"] == mname)
    if(length(infoIndex) == 1 && !is.na(markerinfo[infoIndex, "IlmnStrand"])){
      if(markerinfo[infoIndex, "IlmnStrand"] == "BOT" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        kabadinerGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else if(markerinfo[infoIndex, "IlmnStrand"] == "TOP" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        kabadinerGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else{
        kabadinerGeno[mname, ] <- genos
      }
    }else{
      kabadinerGeno[mname, ] <- genos
    }
  }
}
write.table(kabadinerGeno, "Kabadiner/analysis/genotypes_snp_kabadiner-2.txt", sep="\t", quote=FALSE)


# Our data (Kabadiner horses) 1

horsedata   <- read.table("Kabadiner/input/kabadiner-1.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=2)
genotypes   <- horsedata[15:nrow(horsedata), ]                                      # Row 14 till the end contains genotype data
rownames(genotypes)   <- gsub("-", "_", rownames(horsedata)[15:nrow(horsedata)])
genotypes <- genotypes[which(genotypes[,"GenTrain.Score"] >= 0.6), ]
genotypes <- genotypes[,grep("Top.Alleles", colnames(genotypes))]
colnames(genotypes) <- gsub(".Top.Alleles", "", colnames(genotypes))

kabadinerGeno <- matrix(NA, length(markers), ncol(genotypes), dimnames = list(markers, colnames(genotypes)))

for(x in 1:nrow(kabadinerGeno)){
  mname <- rownames(kabadinerGeno)[x]
  genoIndex <- which(rownames(genotypes) == mname)
  cat(mname, genoIndex, "\n")
  if(length(genoIndex) == 1){
    genos <- unlist(genotypes[genoIndex, ])
    infoIndex <- which(markerinfo[,"Name"] == mname)
    if(length(infoIndex) == 1 && !is.na(markerinfo[infoIndex, "IlmnStrand"])){
      if(markerinfo[infoIndex, "IlmnStrand"] == "BOT" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        kabadinerGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else if(markerinfo[infoIndex, "IlmnStrand"] == "TOP" && markerinfo[infoIndex, "SourceStrand"] == "BOT"){
        kabadinerGeno[mname, ] <- unlist(lapply(lapply(strsplit(genos,""), opposite), paste0, collapse=""))
      }else{
        kabadinerGeno[mname, ] <- genos
      }
    }else{
      kabadinerGeno[mname, ] <- genos
    }
  }
}
write.table(kabadinerGeno, "Kabadiner/analysis/genotypes_snp_kabadiner-1.txt", sep="\t", quote=FALSE)
