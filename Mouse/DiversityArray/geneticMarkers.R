# markersInLab.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Analysis of genetic marker data from Angelika Ackermann

setwd("E:/Mouse/DNA/DiversityArray/")

snpsChr3 <- c("rs29616484", "rs33497409", "rs29591422", "rs48063146", "rs3152105", "rs3694769", "rs29986155", "rs3151486", "rs29993050", "rs3151465", "rs50516118",
              "rs50895669", "rs13477092", "rs30001450", "rs30731601", "rs36659747", "rs31664218", "rs32954388", "rs30106246", "rs3151486", "rs505161118",
              "rs29657774", "rs3151604", "rs3685081", "rs16799508")
snpsChr7 <- c("rs32471611", "rs31371879", "rs32516085", "rs31060001", "rs33140563", "rs33153126", "rs3697227", "rs31373883", "rs32111494", "rs16807448", "rs33182022", 
              "rs47135628", "rs33183784", "rs33185539", "rs31253258", "rs31148458", "rs31308196", "rs51716405", "rs46073185", "rs47605037", "rs33257139", "rs32111494", 
              "rs16807448", "rs33182022","rs47135628")

microSatellites <- rbind(
c("D3Mit46",  "11/22",  3,  30013439,  30013604),
c("D3Mit272", "22/11",  3,  33007369,  33007465),
c("D3Mit21",  "11/22",  3,  37125593,  37125823),
c("D3Mit296", "11/22",  3,  48240972,  48241090),
c("D3Mit183", "22/11",  3,  52176672,  52176810),
c("D3Mit22",  "11/22",  3,  69718176,  69718412),
c("D3Mit312", "22/11",  3, 101400619, 101400741),
c("D3Mit89",  "11/22",  3, 156837105, 156837325))

setwd("E:/Mouse/DNA/")

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

setwd("E:/Mouse/DNA/DiversityArray/")

library(biomaRt)                                              # Biomart
snp.db <- useMart("snp", dataset="mmusculus_snp")             # For mouse SNPs
res.biomart <- getBM(c("refsnp_id", "allele", "chr_name", "chrom_start", "chrom_start", "chrom_start"), filters="snp_filter", values=c(snpsChr3, snpsChr7), mart=snp.db)

res.biomart <- apply(res.biomart, 2, as.character)

colnames(res.biomart)     <- c("markerID", "Allele", "Chr", "Start", "End", "Location")
colnames(microSatellites) <- c("markerID", "Allele", "Chr", "Start", "End", "Location")

markersInLab <- as.data.frame(rbind(microSatellites, res.biomart))
write.table(markersInLab, "Annotation/GeneticMarkers.txt", sep="\t", row.names=FALSE)

