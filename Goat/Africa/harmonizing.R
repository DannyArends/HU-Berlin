setwd("D:/Edrive/Goat/")

data.adaptmap <- read.table(file = "Adaptmap.matrix.txt", sep = "\t")
data.algeria <- read.table(file = "Algeria.matrix.txt", sep = "\t")
data.botswana <- read.table(file = "Botswana.matrix.txt", sep = "\t")
data.ethiopia <- read.table(file = "Ethiopia.matrix.txt", sep = "\t")
data.sudan <- read.table(file = "Sudan.matrix.txt", sep = "\t")
data.swiss <- read.table(file = "Swiss.matrix.txt", sep = "\t")
data.uganda <- read.table(file = "Uganda.matrix.txt", sep = "\t")

snp.annot <- read.table("DNA/SihamRawData/snpinfo.txt", sep = "\t")
snp.annot <- snp.annot[,c("reference", "SNP")]
snp.annot[,"SNP"] <- substring(snp.annot[,"SNP"],2,4)

snps.shared <- c(rownames(data.adaptmap), rownames(data.algeria), rownames(data.botswana), rownames(data.ethiopia), 
                 rownames(data.sudan), rownames(data.swiss), rownames(data.uganda), rownames(snp.annot))
snps.shared <- names(which(table(snps.shared) == 8))

snp.annot <- snp.annot[snps.shared,]

data.adaptmap <- data.adaptmap[snps.shared, ]
data.algeria <- data.algeria[snps.shared, ]
data.botswana <- data.botswana[snps.shared, ]
data.ethiopia <- data.ethiopia[snps.shared, ]
data.sudan <- data.sudan[snps.shared, ]
data.swiss <- data.swiss[snps.shared, ]
data.uganda <- data.uganda[snps.shared, ]

samples.names <- c(colnames(data.adaptmap), colnames(data.algeria), colnames(data.botswana), colnames(data.ethiopia), 
                   colnames(data.sudan), colnames(data.swiss), colnames(data.uganda))
samples.loc <- c(rep("adaptmap", ncol(data.adaptmap)),
                 rep("algeria", ncol(data.algeria)),
                 rep("botswana", ncol(data.botswana)),
                 rep("ethiopia", ncol(data.ethiopia)),
                 rep("sudan", ncol(data.sudan)),
                 rep("swiss", ncol(data.swiss)),
                 rep("uganda", ncol(data.uganda)))
annot.samples <- cbind(samples.names, samples.loc)

data.merged <- cbind(data.adaptmap, data.algeria, data.botswana, data.ethiopia, data.sudan, data.swiss, data.uganda)

opposite <- function(x){
  if(x == "A") return("T")
  if(x == "T") return("A")
  if(x == "C") return("G")
  if(x == "G") return("C")
  stop("should not be here")
}

snp.annot <- cbind(snp.annot, hasRef = apply(snp.annot,1,function(x){
  x[1] %in% unlist(strsplit(as.character(x[2]), "/"))
}))

corrected <- c()
for(x in 1:nrow(snp.annot)){
  alleles <- unlist(strsplit(as.character(snp.annot[x,2]), "/"))
  if(!snp.annot[x,3]) { # Reference is not found in the alleles (flip the SNP)
    alleles[1] <- opposite(alleles[1])
    alleles[2] <- opposite(alleles[2])
    if(alleles[1] == as.character(snp.annot[x, "reference"])){
      corrected <- c(corrected, paste0(alleles[1], "/", alleles[2]))
    }else{
      corrected <- c(corrected, paste0(alleles[2], "/", alleles[1]))
    }
  }else{
    if(alleles[1] == as.character(snp.annot[x, "reference"])){
      corrected <- c(corrected, paste0(alleles[1], "/", alleles[2]))
    }else{
      corrected <- c(corrected, paste0(alleles[2], "/", alleles[1]))
    }
  }
}
snp.annot <- cbind(snp.annot, alleles = corrected)

corrected.matrix <- matrix(NA, nrow(data.merged), ncol(data.merged), dimnames = list(rownames(data.merged), colnames(data.merged)))
for(x in 1:nrow(data.merged)){
  datarow <- as.character(unlist(data.merged[x,]))
  corAlleles <- unlist(strsplit(as.character(snp.annot[x, "alleles"]), "/"))

  gts <- rep(NA, length(datarow))
  for(y in 1:length(datarow)){
    if(is.na(datarow[y])){
      gts[y] <- NA
    }else{
      alleles <- unlist(strsplit(datarow[y], ""))
      if(!all(alleles %in% corAlleles)){
        alleles[1] <- opposite(alleles[1])
        alleles[2] <- opposite(alleles[2])
      }
      if(alleles[1] == alleles[2]) { 
        gts[y] <- paste0(alleles[1], alleles[2])
      }else{
        if(alleles[1] == as.character(snp.annot[x, "reference"])){
          gts[y] <- paste0(alleles[1], alleles[2])
        }else{
          gts[y] <- paste0(alleles[2], alleles[1])
        }
      }
    }
  }
  cat(x, "\n")
  corrected.matrix[x,] <- gts
}

write.table(corrected.matrix, file = "Africa.Swiss.Adaptmap.merged.matrix.txt", sep = "\t", quote=FALSE)
write.table(snp.annot, file = "Africa.Swiss.Adaptmap.snp.annotation.txt", sep = "\t", quote=FALSE)
write.table(annot.samples, file = "Africa.Swiss.Adaptmap.samples.annotation.txt", sep = "\t", quote=FALSE, row.names=FALSE)
