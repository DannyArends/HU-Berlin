setwd("D:/Edrive/Mouse/Texas_Pb")
gtm <- read.table("genotypes_all.txt", sep = "\t", colClasses="character")
map <- read.table("reblasted_map.txt",sep="\t", header=TRUE, row.names=1)
gtm <- gtm[which(rownames(gtm) %in% rownames(map)),]
map <- map[rownames(gtm),]

dim(gtm)
cc11 <- grep("CC011", colnames(gtm))
cc17 <- grep("CC017", colnames(gtm))

cc11.tbl <- apply(gtm[,cc11],1,table)
cc17.tbl <- apply(gtm[,cc17],1,table)

gtm <- gtm[which(unlist(lapply(cc11.tbl, length)) >= 1 & unlist(lapply(cc17.tbl, length)) >= 1),]
map <- map[rownames(gtm),]
table(map[,1])
dim(gtm)

cc11.gts <- unlist(apply(gtm[,cc11],1,function(x){names(which.max(table(x)))}))
cc17.gts <- unlist(apply(gtm[,cc17],1,function(x){names(which.max(table(x)))}))

gtm <- gtm[names(which(cc11.gts != cc17.gts)),]
map <- map[rownames(gtm),]
table(map[,1])
dim(gtm)

missingPmarker <- (apply(apply(gtm, 1, is.na), 2, sum) / ncol(gtm)) * 100
gtm <- gtm[-which(missingPmarker > 10),]
map <- map[rownames(gtm),]
table(map[,1])
dim(gtm)

missingPind <- (apply(apply(gtm, 2, is.na), 2, sum) / nrow(gtm)) * 100
gtm <- gtm[, -which(missingPind > 10)]
map <- map[rownames(gtm),]
table(map[,1])
dim(gtm)
write.table(gtm, "genotypes_all_filtered.txt", sep = "\t", quote = FALSE)

# Start again from here
gtm <- read.table("genotypes_all_filtered.txt", sep = "\t", colClasses="character")
gtm <- gtm[,grep("F2", colnames(gtm))]
map <- map[rownames(gtm),]

chroms <- c(1:19, "X", "Y", "MT")
ordering <- sort(map[,2], index.return=TRUE)
map <- map[ordering$ix,]

mapG <- c()
for(chr in chroms){
  mapG <- rbind(mapG, map[which(map[,1] == chr),])
}
map <- mapG

gtm <- gtm[rownames(mapG),]

seggregates <- names(which(lapply(apply(gtm,1,table), length) > 1))

gtm <- gtm[seggregates,]
map <- map[rownames(gtm), ]
table(map[,1])
dim(gtm)

# Change genotypes to follow normal alphabetical order (since we have no phasing information):
gtm[gtm == "TA"] <- "AT"
gtm[gtm == "GC"] <- "CG"
gtm[gtm == "TC"] <- "CT"
gtm[gtm == "TG"] <- "GT"

table(unlist(gtm))

smallGTgroup <- names(which(unlist(lapply(apply(gtm,1,table),min)) < 20))
for(x in smallGTgroup){
  gt <- names(which.min(table(unlist(gtm[x,]))))
  gtm[x, which(gtm[x,] == gt)] <- NA
  if(min(table(unlist(gtm[x,]))) < 20){
    gt <- names(which.min(table(unlist(gtm[x,]))))
    gtm[x, which(gtm[x,] == gt)] <- NA
  }
}

seggregates <- names(which(lapply(apply(gtm,1,table), length) > 1))
gtm <- gtm[seggregates,]
map <- map[rownames(gtm), ]
table(map[,1])
dim(gtm)

write.table(gtm, "genotypes_F2_filtered_ordered.txt", sep = "\t", quote = FALSE)
write.table(map, "map_ordered.txt", sep = "\t", quote = FALSE)

gts.num <- matrix(NA, nrow(gtm), ncol(gtm), dimnames=list(rownames(gtm), colnames(gtm)))

for(x in 1:nrow(gtm)){
  gttbl <- table(unlist(strsplit(unlist(gtm[x,]), "")))
  minorA <- names(which.min(gttbl)); majorA <- names(which.max(gttbl))
  if(minorA == majorA){ minorA <- names(gttbl)[1];  minorA <- names(gttbl)[2] } 
  hMi <- paste0(minorA,minorA)
  het <- paste0(sort(c(minorA,majorA)), collapse="")
  hMa <- paste0(majorA,majorA)
  gts.num[x, which(gtm[x,] == hMi)] <- -1
  gts.num[x, which(gtm[x,] == het)] <- 0
  gts.num[x, which(gtm[x,] == hMa)] <- 1
}

write.table(gts.num, "genotypes_F2_filtered_ordered_numeric.txt", sep = "\t", quote = FALSE)

# Find some more markers which are weird using correlation between markers
gts.num <- read.table("genotypes_F2_filtered_ordered_numeric.txt", sep = "\t")

gts.corM <- abs(cor(t(gts.num), use="pair"))

for(chr in chroms){
  ii <- which(map[,1] == chr)
  png(paste0("images/chr", chr, ".png"), width=1024, height=1024)
  image(1:length(ii),1:length(ii), gts.corM[ii,ii], main=paste0("Chromosome", chr))
  dev.off()
}

avg.cors <- c()
for(x in 1:nrow(gts.corM)){
  chr <- map[x, "chromsome"]
  onChr <- which(map[,"chromsome"] == chr)
  s <- max(x-5, min(onChr))
  e <- min(x+5, max(onChr))
  avg.cors <- rbind(avg.cors, c(chr, x - min(onChr), s, e,mean(gts.corM[x, s:e])))
}
plot(as.numeric(avg.cors[,5]), col = as.numeric(as.factor(map[, "chromsome"])))

resCor <- cbind(map[which(as.numeric(avg.cors[,5]) < 0.8),],avg.cors[which(as.numeric(avg.cors[,5]) < 0.8),])

badC <- which(rownames(gts.corM) %in% rownames(resCor))
gts.corM <- gts.corM[-badC, -badC]
map.corM <- map[-badC,]
for(chr in chroms){
  ii <- which(map.corM[,1] == chr)
  png(paste0("images/chr", chr, "_Clean.png"), width=1024, height=1024)
  image(1:length(ii),1:length(ii), gts.corM[ii,ii], main=paste0("Chromosome", chr))
  dev.off()
}

