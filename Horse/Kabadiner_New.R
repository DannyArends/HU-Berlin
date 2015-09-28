# Analysis of Kabadiner horse SNP and performance data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015
setwd("E:/Horse/DNA/Kabadiner/")

horsedata1   <- read.table("input/Kabadiner-1.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=2)
horsedata2   <- read.table("input/Kabadiner-2.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=2)
chrInfo      <- read.table("info/chrinfo.txt", sep="\t", header=TRUE)
horsedata    <- cbind(horsedata1, horsedata2)


phenotypes  <- horsedata[1:14, ]                                                    # Row 1 till 13 contains phenotype data
genotypes   <- horsedata[15:nrow(horsedata), ]                                      # Row 14 till the end contains genotype data
genotypes   <- genotypes[which(as.numeric(genotypes[,"GenTrain.Score"]) > 0.6),]    # Keep only high quality calls
map         <- genotypes[,c("Chr","Position")]                                      # Extract the map and split the phenotypes and genotypes

calledgeno  <- genotypes[,grep("GType", colnames(genotypes))]
genotypes   <- genotypes[,c(grep("Top.Alleles", colnames(genotypes)),grep("TOP", colnames(genotypes)))]
phenotypes  <- phenotypes[,grep("GType", colnames(phenotypes))]

colnames(phenotypes) <- gsub(".GType", "", colnames(phenotypes))
colnames(calledgeno) <- gsub(".GType", "", colnames(calledgeno))
colnames(genotypes)  <- gsub(".Top.Alleles", "", colnames(genotypes))
colnames(genotypes)  <- gsub(".TOP", "", colnames(genotypes))

### Data QC
tables          <- apply(genotypes, 1, table)
nonInformative  <- which(unlist(lapply(tables,length)) < 2)                         # Non informative, since we only have 1 genotype
genotypes       <- genotypes[-nonInformative, ]
calledgeno      <- calledgeno[-nonInformative, ]
map             <- map[-nonInformative, ]

notDuplicated <- which(!duplicated(calledgeno))                                     # Duplicated markers
genotypes     <- genotypes[notDuplicated, ]
map           <- map[notDuplicated, ]
cat("Left with", nrow(genotypes), "markers\n")

write.table(map, file="input/cleaned_map.txt", sep = "\t")                          # Save the clean map to disk
write.table(genotypes, file="input/cleaned_genotypes.txt", sep = "\t")              # Save the clean genotypes to disk
write.table(phenotypes, file="input/cleaned_phenotypes.txt", sep = "\t")            # Save the clean phenotypes to disk

## Some basic QG plots of all the individuals relatedness
numgeno <- apply(genotypes, 2, function(x){as.numeric(as.factor(x))})

#numgeno <- numgeno[,-which(colnames(numgeno) %in% c("P3691", "P3611", "P3625","P3525","P3422","P3613","P3423"))]
dendrogram <- as.dendrogram(hclust(dist(t(numgeno), method = "manhattan")))         #TODO: perhaps add the reference horse

cols <- c("red","blue","black","green","orange","gold","pink")
names(cols) <- unique(as.character(phenotypes[3,]))

labelCol <- function(x) {
  if (is.leaf(x)) {
    hclass <- phenotypes[3, attr(x, "label")]                                       # Fetch the class label
    hcol <- cols[hclass]                                                            # Label color
    cat(attr(x, "label"), hclass, hcol, "\n")
    attr(x, "nodePar") <- list(lab.col=hcol)
  }
  return(x)
}
dendrogram.col <- dendrapply(dendrogram, labelCol)
#png("Dendrogram_kabadiner.png", width=1024, height=768)
plot(dendrogram.col)
#dev.off()

# Identical : genotypes[,c("P3623","P3625")] and genotypes[,c("P3422","P3425")]

# Siblings? : genotypes[,c("P3611","P3613")], No too close (only 2000 snps different)
equal    <- sum(apply(genotypes[,c("P3611", "P3613")],1,function(x){ x[1] == x[2] }),na.rm=TRUE)
unequal  <- sum(apply(genotypes[,c("P3611", "P3613")],1,function(x){ x[1] != x[2] }),na.rm=TRUE)

# Should be siblings: genotypes[,c("P3420", "P3611")], YES, they are siblings (~ 1/2 of the SNPs is different)
equal    <- sum(apply(genotypes[,c("P3420", "P3611")],1,function(x){ x[1] == x[2] }),na.rm=TRUE)
unequal  <- sum(apply(genotypes[,c("P3420", "P3611")],1,function(x){ x[1] != x[2] }),na.rm=TRUE)

### Only take the kabadiner horses for QTL mapping
#kabadiner <- colnames(phenotypes)[which(phenotypes[6,] == "Kab")]

#genotypes   <- genotypes[,kabadiner]
#calledgeno  <- calledgeno[,kabadiner]
#phenotypes  <- phenotypes[,kabadiner]

### Chromosome plot to show the location of SNPs
chromosomes  <- as.character(c(1:31, "X", "Y", "MT"))

plot(y=c(0, max(chrInfo[,2])), x=c(1,nrow(chrInfo)), t='n', main="", yaxt="n", xlab="Chromosome", ylab="Length (Mb)", xaxt="n")
abline(h=seq(0, max(chrInfo[,2]), 10000000), col = "lightgray", lty = "dotted")

cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1,lwd=3)
  cnt <<- cnt + 1
})

aa <- apply(map,1,function(x){
  xloc <- match(as.character(x["Chr"]), chromosomes); yloc <- as.numeric(x["Position"])
  points(x=xloc, y=yloc, pch="-",cex=1.0)
})

axis(1, chromosomes, at=c(1:nrow(chrInfo)), las=1)
axis(2, seq(0, max(chrInfo[,2]), 10000000)/1000000, at=seq(0, max(chrInfo[,2]), 10000000), cex.axis=0.7)

Bay     <- paste0("P", c(3509, 3683, 3454, 3674, 3681, 3714, 3632, 3638, 3696, 3728))
DarkBay <- paste0("P", c(3417, 3479, 3640, 3715, 3667))
SealBay <- paste0("P", c(3421, 3571, 3597, 3613, 3007, 3615))

# TODO find a SNP which is different in one group compared to the other 2 groups
inBay <- apply(genotypes[,Bay],1,function(x){length(table(x)) == 1})
inDarkBay <- apply(genotypes[,DarkBay],1,function(x){length(table(x)) == 1})
inSealBay <- apply(genotypes[,SealBay],1,function(x){length(table(x)) == 1})

subsetG <- genotypes[which(inBay & inDarkBay & inSealBay), ]

difBayDarkBay  <- which(apply(subsetG[, c(Bay, DarkBay)],1,function(x){ length(table(x)) == 2 }))
difBaySealBay  <- which(apply(subsetG[, c(Bay, SealBay)],1,function(x){ length(table(x)) == 2 }))
DarkBaySealBay <- which(apply(subsetG[, c(DarkBay, SealBay)],1,function(x){ length(table(x)) == 2 }))

## TODO same analysis now allowing for heterozygous animals in a single group so AA,AB versus BB

## QTL mapping analysis
sex        <- as.factor(unlist(phenotypes["Sex",]))                                     # Has to be corrected for
breed      <- as.factor(unlist(phenotypes["Breed",]))                                   # Has to be corrected for (3)

owner      <- as.factor(unlist(phenotypes["Owner",]))                                   # Might but a lot of levels (41)
birthyear  <- as.factor(unlist(phenotypes["Year",]))                                    # Might but a lot of levels (11)
father     <- as.factor(unlist(phenotypes["Father",]))                                  # Might but a lot of levels (20)
mtDNA      <- as.factor(unlist(phenotypes["mt-Haplo",]))                                # Might but a lot of levels (19)
mtDNA2     <- as.factor(unlist(lapply(strsplit(as.character(mtDNA),""),"[",1)))         # Might (7)

for(pheno in c("Performance")){ #c("Performance","WH","BU","RU")){
  phenotype  <- as.numeric(phenotypes[pheno, ])
  # Sex and mtDNA (class) seem to influence our phenotype
  anova(lm(phenotype ~ sex + breed))
  anova(lm(phenotype ~ sex + mtDNA))

  pvalues <- NULL
  for(x in 1:nrow(genotypes)){
    tryCatch(res <- anova(lm(phenotype ~ sex + breed + as.factor(unlist(genotypes[x,]))))[[5]], error = function(e){ res <<- rep(NA, 3) })
    pvalues <- rbind(pvalues, res[-length(res)])
  }
  colnames(pvalues) <- c("Sex", "Breed", "Marker")
  write.table(cbind(map, -log10(pvalues)), paste0("analysis/all_",pheno,".txt"), sep="\t", row.names=FALSE)
}

# Extract data from the QTL mapping analysis (LOD score above 6, which marker is it)
above <- which(-log10(pvalues[,3]) > 6)
p.adjust(pvalues[,3],"BH")[above]

region.start <- min(map[(which(-log10(pvalues[,3]) > 6)-5):(which(-log10(pvalues[,3]) > 6)+5),"Position"])
region.end <- max(map[(which(-log10(pvalues[,3]) > 6)-5):(which(-log10(pvalues[,3]) > 6)+5),"Position"])

res.bio <- getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene", "external_gene_name","description", "chromosome_name", "start_position", "end_position"), 
filters = c("chromosomal_region", "biotype"), values = list("14:41000000:43000000","protein_coding"), mart = bio.mart)

write.table(res.bio, file="analysis/GenesNearTopSNP.txt",sep="\t",row.names=FALSE)

## Colors for the Manhattan plot
cols <- rep(c("red","blue"),16)
names(cols) <- unique(map[,"Chr"])

## Manhattan plot
for(pheno in c("Performance","WH","BU","RU")){
  pvalues <- read.table(paste0("analysis/",pheno,".txt"), sep="\t",header=TRUE)[,c(3:4)]
  png("Performance_Kabadiner.png", width=1024,height=768)

  plot(pvalues[,2], t = 'h', col=cols[map[,"Chr"]], main="Endurance in Kabadina horses", ylab="LOD score", xlab="SNP")                                                                             # Uncorrected for multiple testing
  abline(h=-log10(0.05/nrow(pvalues)), col="green", lwd=1,lty=2) ; abline(h=-log10(0.01/nrow(pvalues)), col="green", lwd=1,lty=2)      # Significance thresholds
  abline(h=-log10(0.1/nrow(pvalues)), col="orange", lwd=1,lty=2) ; abline(h=-log10(0.01/nrow(pvalues)), col="orange", lwd=1,lty=2)      # Significance thresholds
  dev.off()

#  plot(-log10(p.adjust(pvalues[,2],"BH")),t='h', col=cols[map[,"Chr"]])                                                                 # BH corrected for multiple testing
#  abline(h=-log10(0.05), col="orange", lwd=1,lty=2) ; abline(h=-log10(0.01), col="green", lwd=1,lty=2)                                  # Significance thresholds
}

# Look into our top marker
plot(phenotype ~ as.factor(unlist(genotypes[which(-log10(pvalues[,2]) > 5),])))