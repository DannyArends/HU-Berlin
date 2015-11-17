
input.dir <- "/home/arends/NAS/Cattle/FUGATO_Genotrack_2008-04-01_SNP50K_Chip/Originaldaten/CAU_Kiel_Rawdata/RAWDATA/"
input.archive <- "/home/arends/NAS/Cattle/FUGATO_Genotrack_2008-04-01_SNP50K_Chip/Originaldaten/CAU_Kiel_Rawdata/Archive/"
out.dir <- "/home/arends/NAS/Cattle/FUGATO_Genotrack_2008-04-01_SNP50K_Chip/Originaldaten/CAU_Kiel_Rawdata/analysis/"
files.all <- dir(input.dir)
files.matrix  <- files.all[grepl("_FinalReport_Matrix.txt$", files.all)]
files.snpinfo <- files.all[grepl("_snp.info$", files.all)]

# Read, load and save the genotypes
if(!file.exists(paste0(out.dir,"genotypes.txt"))){
  alldata <- NULL
  for(fname in files.matrix) {
    idx <- which(grepl("ARS-BFGL-BAC-10172", readLines(paste0(input.dir, fname), n=35)))
    if(idx > 30) idx <- (idx-4)
    if(idx < 20) idx <- (idx-2)
    one.data <- read.table(paste0(input.dir, fname), sep = "\t", skip = idx, header = TRUE, check.names = FALSE, row.names = 1, na.strings=c("NA", "--", ""))
    if(is.null(alldata)) {
      alldata <- one.data
    } else {
      alldata <- cbind(alldata, one.data)
    }
    cat("Loaded ", fname, ", nrows: ", nrow(one.data), ", cols: ", ncol(one.data), "\n", sep = "")
  }
  write.table(alldata, file = paste0(out.dir,"genotypes.txt"), quote = FALSE, sep = "\t")
}else{
  alldata <- read.csv(file = paste0(out.dir,"genotypes.txt"), sep = "\t", row.names = 1, header = TRUE, check.names = FALSE)
}

# Read, load and save the snp information
snpinfo <- read.table(paste0(input.dir, "SNP_Map.txt"), sep = "\t", header = TRUE)
rownames(snpinfo) <- snpinfo[,"Name"]
snpinfo <- snpinfo[,-c(1,2,5,7,8,9)]
write.table(snpinfo, paste0(out.dir,"snpinfo.txt"), sep = "\t", quote = FALSE)

# Read, load, update and save the sample information
sampleinfo <- read.csv(paste0(input.dir, "Referenz_Platten1_bis_38.txt"), sep="\t", header = TRUE, row.names=1, colClasses="character")
sampleinfo39 <- read.csv(paste0(input.dir, "reference_p39.txt"), sep="\t", header = TRUE, row.names=1, colClasses="character")
sampleinfo39 <- sampleinfo39[,-8]
sampleinfo <- rbind(sampleinfo, sampleinfo39)
sampleinfo <- cbind(sampleinfo, EU_LOM = NA, NAT_LOM = NA, Origin = NA)

samplearchive <- read.csv(paste0(input.archive, "Reference_update.txt.bak20090709"), sep="\t", header = TRUE, row.names=1, colClasses="character", na.strings=c("NA", "NULL", ""))

for(name in rownames(sampleinfo)){
  idx <- which(rownames(samplearchive) == name)
  if(length(idx) == 1){
    if(is.na(sampleinfo[name, "LOM"])) sampleinfo[name, "LOM"] <- samplearchive[idx, "EU_LOM"]
    sampleinfo[name, "EU_LOM"] <- samplearchive[idx, "EU_LOM"]
    sampleinfo[name, "NAT_LOM"] <- samplearchive[idx, "NAT_LOM"]
    sampleinfo[name, "Origin"] <- samplearchive[idx, "ORIGIN"]
  }
}
write.table(sampleinfo, paste0(out.dir,"sampleinfo.txt"), sep = "\t", quote = FALSE)

length(which(!colnames(alldata) %in% rownames(sampleinfo)))
length(which(!rownames(alldata) %in% rownames(snpinfo)))

# Read, load and save the phenotypes
phenotypes <- read.csv(paste0(input.archive, "Phenotypes_HF.txt.bak091104"), sep = "\t", header = TRUE, na.strings=c("\\N", "--", ""), colClasses="character", check.names = FALSE)

length(which(phenotypes[,"V_RINDID"] %in% sampleinfo[,"LOM"]))    # 1823 genotyped samples are fathers
length(which(phenotypes[,"M_RINDID"] %in% sampleinfo[,"LOM"]))    #    9 genotyped samples are mothers
length(which(phenotypes[,"RINDID"] %in% sampleinfo[,"LOM"]))      # 2432 genotyped samples are animals


# Create numberic data, so we can do a single big dendrogram
numdata <- apply(alldata, 1, function(x){as.numeric(as.factor(x))})
colnames(numdata) <- rownames(alldata)
rownames(numdata) <- colnames(alldata)
distances <- dist(numdata[, sample(ncol(numdata), 1000)], method="manhattan")
clusters <- hclust(distances)
dendro <- as.dendrogram(clusters)

breeds <- unique(sampleinfo[,"breed"])
cols <- rainbow(length(breeds))
names(cols) <- breeds

colLab <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    color <- cols[sampleinfo[a$label, "breed"]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = as.character(color)) # cols[sampleinfo[a$label, "breed"]]
  }
  n
}

dendroC <- dendrapply(dendro, colLab)

png(paste0(out.dir,"dendrogram_all.png"), res = 144, width = 7 * 4096, height = 2400, pointsize = 10)
  op <- par(cex = 0.5)
  plot(dendroC)
dev.off()
