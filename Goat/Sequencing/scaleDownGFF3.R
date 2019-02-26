
### Update the GFF3 files
setwd("D:/Edrive/Goat/DNA/Sequencing")
gff3 <- read.csv("ref_ASM170441v1_top_level.gff3", sep = "\t", comment.char = "#", header = FALSE, colClasses = "character")

header <- readLines("ref_ASM170441v1_top_level.gff3",n=9)

chrN <- rbind(chr6  = c("ENA|CM004567|CM004567.1", "NC_030813.1"),
              chr8  = c("ENA|CM004569|CM004569.1", "NC_030815.1"),
              chr14 = c("ENA|CM004575|CM004575.1", "NC_030821.1"))

gff3 <- gff3[which(gff3[,1] %in% chrN[,2]),]
gff3 <- gff3[which(gff3[,3] != "region"),]
gff3 <- gff3[which(gff3[,3] != "match"),]

for(x in 1:nrow(chrN)){
  gff3[which(gff3[,1] == chrN[x,2]),1] <- chrN[x,1]
}

# Write the updated GFF3
cat(paste0(header, collapse="\n"), "\n", file = "ugff3.gff3")
write.table(gff3, file="ugff3.gff3", sep="\t", append=TRUE,col.names=FALSE, row.names=FALSE, quote=FALSE)
