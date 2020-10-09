setwd("D:/Edrive/Mouse/Texas_Pb")

gt1 <- read.csv("input/Texas_A&M_Threadgill_MURCOMV02_20191118_FinalReport.txt", header=TRUE, skip = 9, sep="\t", na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
gt2 <- read.csv("input/Texas_A&M_Threadgill_MURCOMV02_20200519_FinalReport.txt", header=TRUE, skip = 9, sep="\t", na.strings=c("","-", "na", "NA", "NaN", "X", "x"))
gt3 <- read.csv("input/Texas_A&M_Threadgill_MURCOMV02_20200526_FinalReport.txt", header=TRUE, skip = 9, sep="\t", na.strings=c("","-", "na", "NA", "NaN", "X", "x"))

gts <- rbind(gt1,gt2,gt3)

markers <- unique(gts[,1])
samples <- unique(gts[,2])

# Filter SNPs based on 0.7, as recommended by illumina:
# - https://www.illumina.com/Documents/products/technotes/technote_gencall_data_analysis_software.pdf

gtm <- matrix(NA, length(markers), length(samples), dimnames=list(markers, samples))

for(x in 1:nrow(gts)){
  if(!is.na(gts[x, "GC.Score"]) && gts[x, "GC.Score"] > 0.7){
    gtm[gts[x, 1],gts[x, 2]] <- paste0(gts[x, "Allele1...Forward"], gts[x, "Allele2...Forward"])
  }
}

write.table(gtm, "genotypes_all.txt", sep = "\t", quote = FALSE)

gtm <- read.table("genotypes_all.txt", sep = "\t", colClasses="character")

missingPmarker <- (apply(apply(gtm, 1, is.na), 2, sum) / ncol(gtm)) * 100
gtm <- gtm[-which(missingPmarker > 5),]

missingPind <- (apply(apply(gtm, 2, is.na), 2, sum) / nrow(gtm)) * 100
gtm <- gtm[, -which(missingPind > 5)]

write.table(gtm, "genotypes_all_filtered.txt", sep = "\t", quote = FALSE)

# Start again from here
gtm <- read.table("genotypes_all_filtered.txt", sep = "\t", colClasses="character")
gtmF2 <- gtm[,grep("F2", colnames(gtm))]

map <- read.table("input/SNP_Map.txt",sep="\t", header=TRUE, row.names=2)
map <- map[rownames(gtmF2),]

map[which(map[,2] == 0), "Chromosome"] <- "MT"
map <- map[,c(2,3)]

chroms <- c(1:19, "X", "Y", "MT")
ordering <- sort(map[,2], index.return=TRUE)
map <- map[ordering$ix,]

mapG <- c()
for(chr in chroms){
  mapG <- rbind(mapG, map[which(map[,1] == chr),])
}

gtmF2 <- gtmF2[rownames(mapG),]

seggregates <- names(which(lapply(apply(gtmF2,1,table), length) > 1))

gtmF2 <- gtmF2[seggregates,]
mapG <- mapG[seggregates, ]

write.table(gtmF2, "genotypes_F2_filtered_ordered.txt", sep = "\t", quote = FALSE)
write.table(mapG, "map_ordered.txt", sep = "\t", quote = FALSE)

