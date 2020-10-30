#Re-Blast the probes
setwd("D:/Edrive/Mouse/Texas_Pb")
probes <- read.csv("input/Companion_Animal_ike_20007095X356944_A2_Mus musculus.csv", skip=7, header=TRUE)
map <- read.table("input/SNP_Map.txt",sep="\t", header=TRUE, row.names=2)
setwd("D:/minimugablast")

cat("", file="minimuga.fasta")
for(x in 1:nrow(probes)){
  cat(">", probes[x, "Name"], "\n", sep="", file="minimuga.fasta", append=TRUE)
  cat(unlist(strsplit(probes[x, "AlleleA_ProbeSeq"], "[", fixed=TRUE))[1], "\n",sep="", file="minimuga.fasta", append=TRUE)
}

#makeblastdb -in mm10.fa -dbtype nucl
#blastn -db mm10.fa -ungapped -out minimuga.txt -query minimuga.fasta -outfmt 6

mdata <- read.csv("minimuga.txt",sep="\t",header=FALSE)
mdata <- mdata[-which(mdata[, "V5"] > 4),]

tbls <- table(mdata[,1])
singleMap <- names(which(tbls == 1))



SNPmap <- mdata[which(mdata[,1] %in% singleMap), c(1, 2, 10)]
SNPmap[,2] <- gsub("chr", "", SNPmap[,2])
SNPmap[which(SNPmap[,2] == "M"),2] <- "MT"
colnames(SNPmap) <- c("name", "chromsome", "position")

setwd("D:/Edrive/Mouse/Texas_Pb")
write.table(SNPmap, file="reblasted_map.txt", sep="\t", quote=FALSE, row.names=FALSE)

## Use the probes as map
SNPmap <- probes[, c("Name", "Chr", "MapInfo")]
colnames(SNPmap) <- c("name", "chromsome", "position")
write.table(SNPmap, file="reblasted_map.txt", sep="\t", quote=FALSE, row.names=FALSE)