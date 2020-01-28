# To numeric genotypes

setwd("D:/Edrive/Goat/")
data.all <- read.table(file = "Africa.Swiss.Adaptmap.merged.matrix.txt", sep = "\t", colClasses = "character")

# To numeric genotype
data.num <- matrix(NA, nrow(data.all), ncol(data.all), dimnames=list(rownames(data.all), colnames(data.all)))
for(r in 1:nrow(data.all)){
  alleles <- unlist(strsplit(snp.annot[r, "alleles"], "/"))
  h0 <- paste0(alleles[1], alleles[1])
  he <- paste0(alleles[1], alleles[2])
  h1 <- paste0(alleles[2], alleles[2])
  data.num[r,] <- as.numeric(factor(data.all[r,], levels = c(h0,he,h1)))
  cat(r, "\n")
}
data.num <- (data.num - 2)
write.table(data.num, file = "Africa.Swiss.Adaptmap.merged.num.matrix.txt", sep = "\t", quote=FALSE)

