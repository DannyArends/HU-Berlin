setwd("D:/Edrive/Goat/")

# Our data
data.siham <- read.table("DNA/SihamAnalysis/filtered_snps.txt", sep = "\t")
colnames(data.siham) <- gsub(".", "", colnames(data.siham), fixed=TRUE)

write.table(data.siham, file = "Sudan.matrix.txt", sep = "\t", quote=FALSE)

# Algeria
data.algeria <- read.table("Algeria/AlgerianGoat.ped", sep = "\t", na.strings=c("", "NA", "0", "-"))
map.algeria <- read.table("Algeria/AlgerianGoat.map", sep = "\t")

mat.algeria <- c()
for(i in seq(8, ncol(data.algeria),2)){
  mat.algeria <- rbind(mat.algeria, paste0(data.algeria[,i], data.algeria[,i+1]))
}
mat.algeria[mat.algeria == "NANA"] <- NA
colnames(mat.algeria) <- data.algeria[,1]
rownames(mat.algeria) <- map.algeria[,2]

write.table(mat.algeria, file = "Algeria.matrix.txt", sep = "\t", quote=FALSE)

# Uganda
data.uganda <- read.table("Uganda/goat_filtered_1_ars.ped", sep = " ", na.strings=c("", "NA", "0", "-"))
map.uganda <- read.table("Uganda/goat_filtered_1_ars.map", sep = "\t")

mat.uganda <- c()
for(i in seq(7, ncol(data.uganda),2)){
  mat.uganda <- rbind(mat.uganda, paste0(data.uganda[,i], data.uganda[,i+1]))
}
mat.uganda[mat.uganda == "NANA"] <- NA
colnames(mat.uganda) <- data.uganda[,2]
rownames(mat.uganda) <- map.uganda[,2]

write.table(mat.uganda, file = "Uganda.matrix.txt", sep = "\t", quote=FALSE)

# Botswana
data.botswana <- read.table("Botswana/raw/tswanagoats.ped", sep = " ", na.strings=c("", "NA", "0", "-"))
map.botswana <- read.table("Botswana/raw/tswanagoats.map", sep = " ")

mat.botswana <- c()
for(i in seq(7, ncol(data.botswana),2)){
  mat.botswana <- rbind(mat.botswana, paste0(data.botswana[,i], data.botswana[,i+1]))
}
mat.botswana[mat.botswana == "NANA"] <- NA
colnames(mat.botswana) <- data.botswana[,2]
rownames(mat.botswana) <- map.botswana[,2]

write.table(mat.botswana, file = "Botswana.matrix.txt", sep = "\t", quote=FALSE)

# Egypt
data.egypt <- read.table("Egypt/plinkBarkiGoat.ped", sep = " ", na.strings=c("", "NA", "0", "-"))
map.egypt <- read.table("Egypt/plinkBarkiGoat.map", sep = "\t")

mat.egypt <- c()
for(i in seq(7, ncol(data.egypt),2)){
  mat.egypt <- rbind(mat.egypt, paste0(data.egypt[,i], data.egypt[,i+1]))
}
mat.egypt[mat.egypt == "NANA"] <- NA
colnames(mat.egypt) <- data.egypt[,2]
rownames(mat.egypt) <- map.egypt[,2]

write.table(mat.egypt, file = "Egypt.matrix.txt", sep = "\t", quote=FALSE)

#Ethiopia&Cameron

data.ethiopia <- read.csv("Ethiopia&Cameron/Genotypes_Ethiopia&cameron.txt", row.names=1, sep = "\t", na.strings=c("", "NA", "0", "-", "--"))
write.table(data.ethiopia[,-c(1,2)], file = "Ethiopia.matrix.txt", sep = "\t", quote=FALSE)

#Swiss
data.swiss <- read.table("Dryad/goat_data2_dryad.tped", sep = " ", na.strings=c("", "NA", "0", "-"))
map.swiss <- read.table("Dryad/goat_data2_dryad.tfam", sep = "\t")

mat.swiss <- c()
for(i in seq(5, ncol(data.swiss),2)){
  mat.swiss <- rbind(mat.swiss, paste0(data.swiss[,i], data.swiss[,i+1]))
}
mat.swiss[mat.swiss == "NANA"] <- NA
colnames(mat.swiss) <- data.swiss[,2]
rownames(mat.swiss) <- map.swiss[,2]

write.table(t(mat.swiss), file = "Swiss.matrix.txt", sep = "\t", quote=FALSE)

# Adaptmap
setwd("D:/Edrive/Goat/")
map.adaptmap <- read.table("Adaptmap/adaptmap.map", sep = "\t")
pedData <- readLines("Adaptmap/adaptmap.ped")

mat.adaptmap <- matrix(NA, length(gsub(" ", "", unlist(strsplit(pedData[1], "\t")))[-c(1:6)]), length(pedData))
samples <- c()
breeds <- c()
for(i in 1:length(pedData)){
  gts <- gsub(" ", "", unlist(strsplit(pedData[i], "\t")))
  samples <- c(samples, gts[2])
  breeds <- c(breeds, gts[1])
  mat.adaptmap[,i] <- gts[-c(1:6)]
  cat(i, "\n")
}
mat.adaptmap[mat.adaptmap == "00"] <- NA
colnames(mat.adaptmap) <- samples
rownames(mat.adaptmap) <- map.adaptmap[,2]
write.table(mat.adaptmap, file = "Adaptmap.matrix.txt", sep = "\t", quote=FALSE)
write.table(cbind(samples, breeds), "Adaptmap/breeds.txt", sep = "\t", row.names=FALSE, quote=FALSE)