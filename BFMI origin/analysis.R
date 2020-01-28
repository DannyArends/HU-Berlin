#
# Origin of the BFMI
#

loadChromosome <- function(chromosome = "MT"){
  cat("Loading chromosome ", chromosome, "\n")
  mdata <- readLines(gzfile(paste0("MGP_BFMI_SNPs_Chr", chromosome,".vcf.gz")))
  mdata <- mdata[-grep("##", mdata)]
  mdata <- gsub("#", "", mdata)

  header <- unlist(strsplit(mdata[1], "\t"))
  nrows <- length(mdata) - 1

  mm <- matrix(unlist(strsplit(mdata[-1], "\t")), nrows, length(header), byrow=TRUE)
  colnames(mm) <- header

  complexSNP <- grep(",", mm[, "ALT"])
  mm <- mm[-complexSNP,]

  mm[,10:ncol(mm)] <- t(apply(mm[,10:ncol(mm)], 1, function(x){
    unlist(lapply(strsplit(x, ":"), "[", 1))
  }))

  mm[mm == "./."] <- NA
  mm[mm == "0/0"] <- -1
  mm[mm == "0/1"] <- 0
  mm[mm == "1/1"] <- 1

  mm <- mm[, -which(colnames(mm) %in% c("ID", "QUAL", "REF", "ALT", "FILTER", "INFO", "FORMAT"))]
  positions <- as.numeric(mm[,"POS"])
  mm <- mm[, -c(1:2)]
  mm <- apply(mm, 2, as.numeric)

  # Get rid of SNPs where some strain do not show SNPs
  nonmissing <- which(apply(apply(mm, 1, is.na), 2, sum) == 0)
  mm <- mm[nonmissing,] # No strains are allowed to be missing
  positions <- positions[nonmissing]
  return(list(snps = mm, positions = positions))
}

analyzeRegions <- function(mb, positions){
  cat("Analysis chromosome ", chromosome, "\n")
  possible <- colnames(mb)
  stretchstart <- 1
  x <- 1
  regions <- c()
  while(x <= nrow(mb)) {
    if(x == nrow(mb)){ # Last region
      if(length(possible) == 1){
        regions <- rbind(regions, c(possible, stretchstart, x, positions[stretchstart], positions[x]))
      }else{
        regions <- rbind(regions, c(NA, stretchstart, x, positions[stretchstart], positions[x]))
      }
    }
    if(length(possible) == 1){
      if (!mb[x, possible]){
        regions <- rbind(regions, c(possible, stretchstart, (x-1), positions[stretchstart], positions[x-1]))
        stretchstart <- x
        possible <- colnames(mb)
      }
    } else {
      newpossible <- c()
      for (i in 1:length(possible)) {
        if (mb[x, possible[i]]) {
          newpossible = c(newpossible, possible[i])
        }
      }
      possible <- newpossible
      if(length(possible) == 0){ # we removed all possible strains
        regions <- rbind(regions, c(NA, stretchstart, (x-1), positions[stretchstart], positions[x-1]))
        stretchstart <- x
        possible <- colnames(mb)
      }
    }
    cat("Marker",x,"left", paste0(possible, sep=" "), "\n")
    x = x + 1;
  }
  colnames(regions) <- c("Strain", "IndexS", "IndexE", "PosS", "PosE")
  return(regions)
}

setwd("D:/Edrive/Mouse/BFMIorigin")
#chromosomes <- c(as.character(1:19), "X", "Y", "MT")
chromosomes <- c("3")
origin <- vector("list", length(chromosomes))
names(origin) <- chromosomes 
for(chromosome in chromosomes){
  res <- loadChromosome(chromosome)
  mm <- res$snps

  clustering <- hclust(dist(t(mm), method = "manhattan"))
  plot(clustering, main = paste0("Chromosome ", chromosome), ylab=paste0("Manhattan distance ",nrow(mm)," SNPs"), xlab = "MGP - mouse strain")

  mb <- mm[,-1] == mm[,1]
  
  regions <- analyzeRegions(mb, res$positions)
  
  write.table(regions, sep="\t", file=paste0("regions_", chromosome, ".txt"), quote=FALSE, row.names=FALSE)
  cat("Done chromosome ", chromosome, "\n")
  gc()
}

regions <- vector("list", length(chromosomes))
names(regions) <- chromosomes 
for(x in 1:length(chromosomes)){
  regions[[chromosomes[x]]] <- read.table(file=paste0("regions_", chromosomes[x], ".txt"), sep="\t", header=TRUE)
}

### Full genome plot
library(RColorBrewer)
colorscheme <- c("black", brewer.pal(12, "Set3"))
maxX <- max(unlist(lapply(regions, function(x){return(max(x[,5]))})))
plot(c(0, maxX), y = c(0, length(chromosomes)), t = 'n', xaxt='n', yaxt='n', xlab="Position (Mb)", ylab="Chromosome")
chr <- 1
for (x in 1:length(regions)) {
  region <- regions[[x]][which(!is.na(regions[[x]][,1])),]
  for(i in 1:nrow(region)){
    color <- colorscheme[as.numeric(factor(as.character(region[i,1]), levels = colnames(mb)))]
    rect(region[i,4], x - 0.3, region[i,5],x + 0.3, col =color, border = NA)
  }
  rect(1, x - 0.3, max(region[,5]),x + 0.3, col = NA)
}
axis(1, at = seq(0, maxX, 10000000), seq(0, maxX, 10000000) / 1000000)
axis(2, at = 1:length(regions), chromosomes,las=2)
legend("right", fill = colorscheme, colnames(mb))

zoomTo <- function(chromosome = 1, start = 1, end = 30000000, addGenes = TRUE){
  library(RColorBrewer)
  colorscheme <- c("black", brewer.pal(12, "Set3"))
  plot(c(start, end), y = c(0, 1), t = 'n', xaxt='n', yaxt='n', xlab="Position (Mb)", ylab="", main=paste0("Chromosome", chromosome, ":", start, "-", end))
  region <- regions[[chromosome]][which(!is.na(regions[[chromosome]][,1])),]
  for(i in 1:nrow(region)){
    color <- colorscheme[as.numeric(factor(as.character(region[i,1]), levels = colnames(mb)))]
    rect(region[i,4], 0.3, region[i,5], 0.5, col = color, border = NA)
  }
  rect(1, 0.3, max(region[,5]), 0.5, col = NA)
  axis(1, at = seq(start, end, (end-start) /10), seq(start, end, (end-start) /10))
  legend("topright", fill = c(colorscheme, "white"), c(colnames(mb), "Undetermined"), ncol=2)
  box()
  if (addGenes) {
    require(biomaRt)
    bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    region <- paste0(chromosome, ":", start, ":", end)
    genesInRegion <- NULL
    res.biomart <- getBM(attributes = c("start_position", "end_position", "strand", "mgi_symbol"), 
                     filters = c("chromosomal_region", "biotype"), values = list(region, "protein_coding"), mart = bio.mart)
    for(x in 1:nrow(res.biomart)){
      strand <- res.biomart[x, 3]
      pos <- 0.6 + (x%%3 / 20)
      if(strand == -1) pos <- 0.2 - (x%%3 / 20)
      rect(res.biomart[x, 1], pos, res.biomart[x, 2], pos + 0.05)
      text((res.biomart[x, 1] + res.biomart[x, 2]) / 2, pos+0.025, res.biomart[x, 4],)
    }
  }
}
zoomTo(19, 4886878, 4906628) # Bbs1
zoomTo(8, 94060000, 94100000) # Bbs2 - NZW_LacJ
zoomTo(9, 59321990, 59353508) # Bbs4
zoomTo(2, 69647171, 69667571) # Bbs5 - AKR/J
zoomTo("3", 36500000, 36700000) # Bbs7 - AKR/J
# Add deletion information


zoomTo(12, 98920574, 98983238) # Bbs8 - SPRET / CAST / NZO
zoomTo(9, 22475715, 22888280) # Bbs9 - C3H_HeJ
zoomTo(10, 111298679, 111301727) # Bbs10
zoomTo(19, 53929861, 53944627) # Bbip1 - C57BL_6NJ
zoomTo(16, 59612949, 59639391) # Arl6 - NZO

zoomTo(8, 92900000, 93500000)
zoomTo("MT", 0, 17000)





library(ape)
phylo <- as.phylo(clustering)
rrphylo <- root(phylo, "CAST_EiJ", resolve.root = TRUE)
plot(rrphylo)

equality <- (mm[, -1] == mm[,"BFMI860-S12"])


wsize <- 10000
binM <- c()
startpos <- min(position)
while (startpos < max(position)) {
  inBIN <- which(position > startpos & position < startpos + wsize)
  if (length(inBIN) > 1) {
    distances <- dist(t(mm[inBIN,]), method = "manhattan", diag = TRUE, upper = TRUE)
    rowdata <- c(startpos, startpos+wsize, length(inBIN), (length(inBIN) * 2) - as.matrix(distances)["BFMI860-S12",])
    plot(hclust(distances))
    Sys.sleep(5)
    binM <- rbind(binM, rowdata)
  }else{
    binM <- rbind(binM, c(startpos, startpos+wsize, 0, rep(NA,14)))
  }
  startpos <- startpos + (wsize / 10)
}
rownames(binM) <- seq(1:nrow(binM))

binP <- apply(binM[, -c(1:3)],1,function(x){
  round(100 * x / sum(x), 2)
})

undetermined <- as.numeric(c(which(apply(apply(binP,2,is.na), 2, sum) > 0), which(apply(apply(binP,2,is.nan), 2, sum) > 0)))
if(length(undetermined) > 0) binP <- binP[,-undetermined]

mids <- apply(binM[as.numeric(colnames(binP)), 1:2],1,sum) /2

plot(c(min(mids),max(mids)), c(0,100), t = 'n')
p <- rep(0, ncol(binP))
for(x in 1:nrow(binP)){
  polygon(x = c(mids,rev(mids)), y = c(p+binP[x,], rev(p)), col = rainbow(13)[x])
  p <- p+binP[x,]
}

# For each bin, determine which strain the BFMI is most similar to
colnames(binM)[3:ncol(binM)][apply(binM[,3:ncol(binM)] / binM[,"BFMI860-S12"],1,which.max)]

# For each bin, plot the similarity between strains
plot(x = c(0, nrow(binM)), y = c(0,1), t='n')
points((binM[,3:ncol(binM)] / binM[,"BFMI860-S12"])[,"C57BL_6NJ"], t = 'l', col="blue")
points((binM[,3:ncol(binM)] / binM[,"BFMI860-S12"])[,"PWK_PhJ"], t = 'l', col="green")
points((binM[,3:ncol(binM)] / binM[,"BFMI860-S12"])[,"CAST_EiJ"], t = 'l', col="red")

# For each bin, plot the similarity between strains
sort(table(colnames(binM)[3:ncol(binM)][apply(binM[,3:ncol(binM)] / binM[,"BFMI860-S12"],1,which.max)]), decreasing = TRUE)
