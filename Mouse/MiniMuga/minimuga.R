setwd("D:/Edrive/Mouse/DNA/MiniMuga")

mdata <- read.csv("miniMUGA-Genotypes.csv", sep=',', na.strings="N", colClasses="character", row.names=1)
mdata <- mdata[mdata[,"Chromosome"] != 0,]
colnames(mdata)[2] <- "Position"
colnames(mdata) <- gsub("BFMI.Gudrun.Brockmann_", "", colnames(mdata))
samples <- rbind(c("M11113", "BFMI860"), c("M28334", "BFMI860-12"), c("M7194", "BFMI861-S1"), c("M7236","BFMI861-S2"))
colnames(mdata)[3:6] <- samples[,2]

chromosomes <- c(1:19, "X", "Y", "PAR", "MT")
chrinfo <- NULL
for(chr in chromosomes){
  chrinfo <- rbind(chrinfo, c(chr, max(as.numeric(mdata[which(mdata[, "Chromosome"] == chr), "Position"]))))
}
colnames(chrinfo) <- c("Chromosome", "Position")
maxy <-  max(as.numeric(chrinfo[,"Position"])) * 1.15

neq <- apply(mdata[,c("BFMI861-S1", "BFMI861-S2")], 1, function(x){
  x[1] != x[2]
})

mdata <- mdata[which(neq), ]

plot(c(1,length(chromosomes)), y = c(0,maxy), t = 'n', xaxt='n', yaxt='n', xlab = "Chromosome", ylab = "Position (Mb)", main="S1 versus S2")
cnt <- 1
apply(chrinfo, 1, function(x){
  cat(as.numeric(x[2]), "\n")
  lines(c(cnt,cnt), c(0,as.numeric(x[2])))
  onChr <- mdata[mdata[,"Chromosome"] == x[1],]
  points(rep(cnt, nrow(onChr)), as.numeric(onChr[,"Position"]), pch="-")
  cnt <<- cnt + 1
})

axis(1, at = 1:length(chromosomes), chromosomes)
axis(2, at = seq(0, maxy, 10000000), seq(0, maxy, 10000000) / 1000000, las=2)

diffs <- diff(as.numeric(mdata[,"Position"]))

hist(diffs[-which(diffs < 0 | diffs > 5000000)])


### microsatelite markers available in the lab
microSatellites <- rbind(
c("D3Mit46",  "11/22",  3,  30013439,  30013604), c("D3Mit272", "22/11",  3,  33007369,  33007465),
c("D3Mit21",  "11/22",  3,  37125593,  37125823), c("D3Mit296", "11/22",  3,  48240972,  48241090),
c("D3Mit183", "22/11",  3,  52176672,  52176810), c("D3Mit22",  "11/22",  3,  69718176,  69718412),
c("D3Mit312", "22/11",  3, 101400619, 101400741), c("D3Mit89",  "11/22",  3, 156837105, 156837325))

setwd("D:/Edrive/Mouse/DNA/")
inLab <- read.table("Annotation/LabGeneticMarkers.txt", sep="\t")
mgiMarkerInfo <- read.table("Annotation/MGImarkers.txt", sep="\t", header=TRUE, colClasses=c("character"))

markers <- as.character(inLab[grep("Mit", inLab[,1]),1])

for(m in markers){
  mgiID <- which(mgiMarkerInfo[,"Symbol"] %in% m)
  if(length(mgiID) > 0){
    if(!(m %in% microSatellites[,1]) && mgiMarkerInfo[mgiID,"Start"] != ""){
      microSatellites <- rbind(microSatellites, c(m, NA, as.character(mgiMarkerInfo[mgiID,c("Chromosome","Start","End")])))
    }
  }
}

microSatellites <- cbind(microSatellites, round((as.numeric(microSatellites[,4]) + as.numeric(microSatellites[,5]))/2, d=0))

colnames(microSatellites) <- c("markerID", "Alleles", "Chr", "Start", "Stop", "Location")

### Differences between S1 and S2
setwd("D:/Edrive/Mouse/DNA/DiversityArray/")
isDuplicated  <- c("BFMI861-S2", "BFMI861-S1")                                                    # Strains measured multiple times

SNPdata <- read.table("Analysis/measurementsAtlas_annotated.txt", sep="\t", header=TRUE)
cat("Starting with", nrow(SNPdata), "valid SNPs from the mouse diversity array\n")

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                           # Load the annotation
validAnnot <- match(colnames(SNPdata)[9:ncol(SNPdata)], rownames(annotation))                     # Annotation that matches the mice we have in our SNP data
annotation <- annotation[validAnnot, ]

inconsistentSNPs <- NULL
amount <- NULL
for(mouseLine in isDuplicated){                                                                   # Find the inconsistent SNPs between the strains
  mouseNames <- rownames(annotation[which(annotation[,"Line"] == mouseLine),])
  unEqual <- which(apply(SNPdata[, mouseNames], 1, function(x){x[1] != x[2]}))
  cat("Found", length(unEqual), "differences between duplicates of", mouseLine,"\n")
  amount <- c(amount, length(unEqual))
  inconsistentSNPs <- c(inconsistentSNPs, unEqual)
}
inconsistentSNPs <- unique(inconsistentSNPs)

percentage <- paste0("(", round(length(inconsistentSNPs)/nrow(SNPdata)*100, d=1), "%)")           # % of inconsistent SNPs
cat("Removing", length(inconsistentSNPs), percentage, "inconsistently genotyped SNPs\n")
SNPinconsistent <- SNPdata[inconsistentSNPs, ]                                                    # Inconsistent SNPs
SNPdata <- SNPdata[-inconsistentSNPs, ]                                                           # Throw away the inconsistent SNPs

###### Comparison: BFMI861-S2 versus BFMI861-S1 ######
bfmiS1 <- rownames(annotation[which(annotation[,"Line"] == "BFMI861-S1"),])[1]
bfmiS2 <- rownames(annotation[which(annotation[,"Line"] == "BFMI861-S2"),])[1]

res <- apply(SNPdata[,c(bfmiS1, bfmiS2)],1,function(x){
  if(is.na(x[1]) || is.na(x[2])) return(FALSE)
  if(x[1] != x[2]) return(TRUE)
  return(FALSE)
})

cat(sum(as.numeric(res),na.rm=TRUE),"Out of",length(res),"\n")
combined <- cbind(SNPdata[which(res),1:8], SNPdata[which(res), bfmiS1], SNPdata[which(res), bfmiS2])      # Create the ouput (subset the whole SNP arrays)
colnames(combined)[9:10] <- c("BFMI861-S1", "BFMI861-S2")

### Code to detect regions
findRegions <- function(snps, snpInRegion = 35, maxDistance = 4000000){
  chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

  regions <- NULL
  for(chr in chromosomes){
    onchr  <- snps[snps[,"Chr"] == chr,]
    if(nrow(onchr) > snpInRegion){
      sorted <- sort(as.numeric(onchr[,"Location"]), index.return=TRUE)
      onchr  <- onchr[sorted$ix,]
      sloc   <- 0
      snpcnt <- 0
      for(snp in 1:(nrow(onchr)-1)){
        cloc <- onchr[snp, "Location"]
        nloc <- onchr[snp+1, "Location"]
        if(as.numeric(nloc) - as.numeric(cloc) < maxDistance){
          if(snpcnt == 0) sloc <- cloc
          snpcnt <- snpcnt + 1
        }else{
          if(snpcnt > snpInRegion){
            regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
            cat("New region:", chr, sloc, cloc, snpcnt,"\n")
          }
          snpcnt <- 0
        }
      }
      if(snpcnt > snpInRegion){
        regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
        cat("New region:",chr, sloc, cloc, snpcnt,"\n")
      }
    }
  }

  colnames(regions) <- c("Chr","Start","End","SNPs")
  return(regions)
}

regions <- findRegions(combined, 50, 4000000)



