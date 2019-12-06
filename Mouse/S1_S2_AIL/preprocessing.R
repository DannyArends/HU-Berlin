
setwd("D:/Edrive/Mouse/S1_S2")

locusxdnaheader <- unlist(strsplit(readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv", n=16)[16],","))
locusxdnasnps <- unlist(strsplit(readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv", n=18)[18],","))

locusxdna <- readLines("Humboldt_Univ_Zu_Berlin_MURGIGV01_20191205_LocusXDNA.csv")[-c(1:22)]
splitted <- strsplit(locusxdna, ",")

calls <- matrix(NA, length(locusxdna) / 2, length(splitted[[1]]))
scores <- matrix(NA, length(locusxdna) / 2, length(splitted[[1]]))
for(x in 1:nrow(mm)){
  if(x %% 2 == 1) calls[x/2,] <- splitted[[x]]
  if(x %% 2 == 0) scores[x/2,] <- splitted[[x]]
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

gts <- t(gts)
idx <- which(apply(gts,1, function(x){sum(is.na(x)) == length(x)}))
gts <- gts[-idx,]

idx <- which(apply(gts,1,function(x){length(table(x)) == 1}))
gts <- gts[-idx,]

gts <- gts[, -ncol(gts)]
write.table(gts, "genotypes.raw.txt", sep="\t", quote=FALSE)

map <- read.table("snp_map.karl.txt", sep = ",", header = TRUE, row.names=1)
map <- map[rownames(gts),]

chrs <- 1:21 
names(chrs) <- c(1:19, "X", "Y")

plot(c(1,21), c(0,200000000), t = 'n', xaxt = "n", las= 2, ylab = "Pos", xlab = "Chr")
aa <- apply(map, 1, function(r) { points(x = chrs[r[1]], y = r[2], pch = "-"); })
axis(1, at = chrs, names(chrs))




