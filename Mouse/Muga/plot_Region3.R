#plot


library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
bm.above <- getBM(attributes = c("ensembl_gene_id","mgi_symbol", "chromosome_name", "start_position", "end_position"), 
filters = c("chromosomal_region", "biotype"), values = list("3:32500000:40000000","protein_coding"), mart = bio.mart)



setwd("E:/Mouse/DNA/MegaMuga/")                                                                                                                                   # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt",             sep="\t", check.names=FALSE, colClasses="character")                                              # Normal A, H, B genotypes
phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)

F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals

genotypes <- genotypes[, F2]
phenotypes <- phenotypes[F2, ]


setwd("E:/Mouse/ClassicalPhenotypes/AIL")
qtl42Fat     <- read.table("Analysis/qtls_fat42.txt",   sep="\t")
qtl42Lean    <- read.table("Analysis/qtls_lean42.txt",   sep="\t")
qtl42FatLean <- read.table("Analysis/qtls_fatDlean42.txt",   sep="\t")

submap <- map
submap <- map[rownames(qtl42Fat[qtl42Fat[,"marker"] > 1,]), ]
submap <- submap[submap[,"Chr"] == 3,]                                    #  & 
submap <- submap[rownames(submap) %in% rownames(qtl42Fat),]

op <- par(mfrow = c(2,1))
plot(c(0, max(as.numeric(submap[,2]))), y = c(0, 65), t = 'n', ylab = "LOD", xlab = "Chromosome 3", xaxt='n', las = 2, main= "MRI measurements (Day 42)")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = qtl42Fat[rownames(submap), "marker"], t = 'l', col = "red")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = qtl42Lean[rownames(submap), "marker"], t = 'l', col = "blue")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = qtl42FatLean[rownames(submap), "marker"], t = 'l', col = "green")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = rep(-1.3, nrow(submap)), pch="|", col = "black", cex=0.5)
legend("topright", c("MRI fat mass","MRI lean mass", "Fat / Lean mass"), col=c("red","blue","green"), lwd=1)
axis(1, at=seq(0, max(as.numeric(submap[,2])), 10000000), seq(0, max(as.numeric(submap[,2])), 10000000) / 1000000)

submap <- submap[as.numeric(submap[,"Mb_NCBI38"]) > 32500000 &  as.numeric(submap[,"Mb_NCBI38"]) < 40000000,]
op <- par(mar = c(5, 4, 0, 2) + 0.1)

plot(c(32500000, 40000000), y = c(-15, 65), t = 'n', ylab = "LOD", xlab = "Chromosome 3: 30mb - 45 mb", xaxt='n', las = 2, main= "")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = qtl42Fat[rownames(submap), "marker"], t = 'l', col = "red")
#points(x = as.numeric(submap[,"Mb_NCBI38"]), y = qtl42Lean[rownames(submap), "marker"], t = 'l', col = "blue")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = qtl42FatLean[rownames(submap), "marker"], t = 'l', col = "green")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = rep(4, nrow(submap)), pch="|", col = "black", cex=0.5)
legend("topright", c("MRI fat mass", "Fat / Lean mass","Gene", "Gm/Riken Gene"), col=c("red","green","gray", "orange"), lwd=c(1,1,4,4))
axis(1, at=seq(32500000, 40000000, 2500000), seq(32500000, 40000000, 2500000) / 1000000)

for(x in 1:nrow(bm.above)){
  offset <- -(4 + x %% 9)
  col <- "gray"
  if(grepl("Gm", bm.above[x,"mgi_symbol"]) || grepl("Rik", bm.above[x,"mgi_symbol"])) col <- "orange"
  lines(x=c(bm.above[x,"start_position"], bm.above[x,"end_position"]), y = c(offset, offset), lwd = 4, col=col)
}


ind <- ncol(genotypes[rownames(submap),])
mar <- nrow(genotypes[rownames(submap),])

sortedByFat <- rownames(phenotypes[sort(phenotypes[,"mri42d_fat"] / phenotypes[,"mri42d_lean"],index.return=TRUE)$ix,])

op <- par(mfrow = c(2,1))
plot(c(30000000, 50000000), y = c(0, 100), t = 'n', ylab = "LOD", xlab = "Region Chr3 20-50 mb", xaxt='n', las = 2, main= "Fat (Day 42)")
points(x = as.numeric(submap[,"Mb_NCBI38"]), y = qtl42Fat[rownames(submap), "marker"], t = 'h')

plot(c(30000000, 50000000), y = c(0, ind), t = 'n', ylab = "Individual", xlab = "Region Chr3 20-50 mb", xaxt='n', las = 2)
sortedgeno <- genotypes[,sortedByFat]


chr3marker <- "UNC5048297"

heavy <- names(which(unlist(sortedgeno[chr3marker,]) == "A"))
other <- names(which(unlist(sortedgeno[chr3marker,]) != "A"))

sortedgeno <- cbind(sortedgeno[,other], sortedgeno[,heavy])
sortedgeno[sortedgeno == "A"] <- 1
sortedgeno[sortedgeno == "H"] <- 2
sortedgeno[sortedgeno == "B"] <- 3

for(i in 1:ind){
  x = as.numeric(submap[, "Mb_NCBI38"])
  mcol <- c("red", "yellow","gray")[as.numeric(as.factor(unlist(sortedgeno[rownames(submap),i])))]
  points(x, rep(i, mar), col = mcol, pch=17, cex=0.5)
}
