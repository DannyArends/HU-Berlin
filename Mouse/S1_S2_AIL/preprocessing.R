setwd("D:/Edrive/Mouse/S1_S2")

locusxdnaheader <- unlist(strsplit(readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv", n=16)[16],","))
locusxdnasnps <- unlist(strsplit(readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv", n=18)[18],","))

locusxdna <- readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv")[-c(1:22)]
splitted <- strsplit(locusxdna, ",")

calls <- matrix(NA, length(locusxdna) / 2, length(splitted[[1]]))
scores <- matrix(NA, length(locusxdna) / 2, length(splitted[[1]]))
for(x in 1:length(splitted)) {
  if(x %% 2 == 1) calls[ceiling(x/2),] <- splitted[[x]]
  if(x %% 2 == 0) scores[ceiling(x/2),] <- splitted[[x]]
}

markers <- locusxdnaheader[4:length(locusxdnaheader)]

colnames(calls) <- c("Label", "plateWell", "Date","oligoPoolId","bundleId", "status", "Type", "Nas", markers)
colnames(scores) <- c("Label", "plateWell", "Date","oligoPoolId","bundleId", "status", "Type", "Nas", markers)

gts <- calls[,markers]
rownames(gts) <- gsub("V 888-", "AIL", calls[, "Label"])
qual <- apply(scores[,markers],2,as.numeric)
rownames(qual) <- gsub("V 888-", "AIL", calls[, "Label"])

gts[qual < 0.7] <- NA
gts[gts == "U"] <- NA

# Write out the raw genotypes
gts <- t(gts)
write.table(gts, "genotypes.raw.txt", sep="\t", quote=FALSE)

gts <- read.table("genotypes.raw.txt", sep="\t")

# Groups with less than 10 individuals are set to missing
tbls <- apply(gts, 1, table)
for(x in 1:length(tbls)){
  for(gt in names(tbls[[x]])){
    if(tbls[[x]][gt] < 10){
      gts[x, which(gts[x,] == gt)] <- NA
    }
  }
}

# All missing
idx <- which(apply(gts,1, function(x){sum(is.na(x)) == length(x)}))
gts <- gts[-idx,]

# Not segregating
idx <- which(apply(gts,1,function(x){length(table(x)) == 1}))
gts <- gts[-idx,]

# More than 10 % missing data
ismissing <- apply(apply(gts, 1, is.na),2,sum)
tooMuchMissing <- names(which((ismissing / ncol(gts)) > 0.1))
gts <- gts[-which(rownames(gts) %in% tooMuchMissing),]

tbls <- apply(gts, 1, table)

map <- read.table("snp_map.karl.txt", sep = ",", header = TRUE, row.names=1)
map <- map[rownames(gts),]

chrs <- 1:21 
names(chrs) <- c(1:19, "X", "Y")

plot(c(1,21), c(0,200000000), t = 'n', xaxt = "n", las= 2, ylab = "Position (mb)", xlab = "Chr", yaxt = 'n')
aa <- apply(map, 1, function(r) { points(x = chrs[r[1]], y = r[2], pch = "-"); })
axis(1, at = chrs, names(chrs))
axis(2, at = seq(0,200000000, 25000000), seq(0,200000000, 25000000)/1000000)

write.table(gts, "genotypes.cleaned.txt", sep="\t", quote=FALSE)
write.table(map, "map.cleaned.txt", sep="\t", quote=FALSE)
