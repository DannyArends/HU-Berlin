# Analysis of Kabadiner horse SNP and performance data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015
setwd("E:/Horse/DNA/Kabadiner/")

horsedata   <- read.table("input/Kabadiner.txt", header=TRUE, sep = "\t",na.strings=c("--", "x", "unknown", ""), colClasses="character", row.names=2)
chrInfo     <- read.table("info/chrinfo.txt", sep="\t", header=TRUE)

phenotypes  <- horsedata[1:13, ]                                                    ## Row 1 till 13 contains phenotype data
genotypes   <- horsedata[14:nrow(horsedata), ]                                      ## Row 14 till the end contains genotype data
genotypes   <- genotypes[which(as.numeric(genotypes[,"GenTrain.Score"]) > 0.6),]    ## Keep only high quality calls

## Extract the map and split the phenotypes and genotypes
map <- genotypes[,c("Chr","Position")]

calledgeno  <- genotypes[,grep("GType", colnames(genotypes))]
genotypes   <- genotypes[,grep("Top.Alleles", colnames(genotypes))]

phenotypes  <- phenotypes[,grep("GType", colnames(phenotypes))]

colnames(phenotypes) <- gsub(".GType", "", colnames(phenotypes))
colnames(calledgeno) <- gsub(".GType", "", colnames(calledgeno))
colnames(genotypes)  <- gsub(".Top.Alleles", "", colnames(genotypes))

## Some basic QG plots of all the individuals relatedness
numgeno <- apply(genotypes, 2, function(x){as.numeric(as.factor(x))})
dendrogram <- as.dendrogram(hclust(dist(t(numgeno), method = "manhattan")))

cols <- c("red","blue","black")
names(cols) <- unique(as.character(phenotypes[6,]))

labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    hclass <- phenotypes[6,label]
    hcol <- cols[hclass]
    cat(label, hclass, hcol, "\n")
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- list(lab.col=hcol)
  }
  return(x)
}
dendrogram.col <- dendrapply(dendrogram, labelCol)

plot(dendrogram.col)

# identical : genotypes[,c("P3623","P3625")]
# identical : genotypes[,c("P3422","P3425")]

# Siblings? : genotypes[,c("P3611","P3613")], NO too close (only 2000 snps different)
equal    <- sum(apply(genotypes[,c("P3611", "P3613")],1,function(x){ x[1] == x[2] }),na.rm=TRUE)
unequal  <- sum(apply(genotypes[,c("P3611", "P3613")],1,function(x){ x[1] != x[2] }),na.rm=TRUE)

# Should be siblings: genotypes[,c("P3420", "P3611")], YES, they are siblings (~ 1/2 of the SNPs is different)
equal    <- sum(apply(genotypes[,c("P3420", "P3611")],1,function(x){ x[1] == x[2] }),na.rm=TRUE)
unequal  <- sum(apply(genotypes[,c("P3420", "P3611")],1,function(x){ x[1] != x[2] }),na.rm=TRUE)

### Only take the kabadiner horses for QTL mapping
kabadiner <- colnames(phenotypes)[which(phenotypes[6,] == "Kab")]

genotypes   <- genotypes[,kabadiner]
calledgeno  <- calledgeno[,kabadiner]
phenotypes  <- phenotypes[,kabadiner]

### Data QC
tables          <- apply(genotypes, 1, table)
nonInformative  <- which(unlist(lapply(tables,length)) < 2)   ## Non informative, since we only have 1 genotype
genotypes       <- genotypes[-nonInformative, ]
calledgeno      <- calledgeno[-nonInformative, ]
map             <- map[-nonInformative, ]

notDuplicated <- which(!duplicated(calledgeno))                ## Duplicated markers
genotypes     <- genotypes[notDuplicated, ]
map           <- map[notDuplicated, ]
cat("Left with", nrow(genotypes), "markers\n")

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

## QTL mapping analysis
sex        <- as.factor(unlist(phenotypes[3,]))                                         # Has to be corrected for
owner      <- as.factor(unlist(phenotypes[2,]))                                         # Might but a lot of levels (41)
birthyear  <- as.factor(unlist(phenotypes[4,]))                                         # Might but a lot of levels (11)
father     <- as.factor(unlist(phenotypes[5,]))                                         # Might but a lot of levels (20)
mtDNA      <- as.factor(unlist(phenotypes[13,]))                                        # Might but a lot of levels (19)
mtDNA2     <- as.factor(unlist(lapply(strsplit(as.character(mtDNA),""),"[",1)))         # Might (7)

for(pheno in c("Performance","WH","BU","RU")){
  phenotype  <- as.numeric(phenotypes[pheno, ])
  # Sex and mtDNA (class) seem to influence our phenotype
  anova(lm(phenotype ~ sex + mtDNA))

  pvalues <- NULL
  for(x in 1:nrow(genotypes)){
    tryCatch(res <- anova(lm(phenotype ~ sex + as.factor(unlist(genotypes[x,]))))[[5]], error = function(e){ res <<- rep(NA, 3) })
    pvalues <- rbind(pvalues, res[-length(res)])
  }
  colnames(pvalues) <- c("sex", "marker")
  write.table(cbind(map, -log10(pvalues)), paste0("analysis/",pheno,".txt"), sep="\t", row.names=FALSE)
}

# Extract data from the QTL mapping analysis (LOD score above 6, which marker is it)
map[which(-log10(pvalues[,2]) > 6),]

## Colors for the Manhattan plot
cols <- rep(c("black","orange"),16)
names(cols) <- unique(map[,"Chr"])

## Manhattan plot
for(pheno in c("Performance","WH","BU","RU")){
  pvalues <- read.table(paste0("analysis/",pheno,".txt"), sep="\t",header=TRUE)[,c(3:4)]
  plot(pvalues[,2], t = 'h', col=cols[map[,"Chr"]])                                                                             # Uncorrected for multiple testing
  abline(h=-log10(0.05/nrow(pvalues)), col="orange", lwd=1,lty=2) ; abline(h=-log10(0.01/nrow(pvalues)), col="green", lwd=1,lty=2)      # Significance thresholds

#  plot(-log10(p.adjust(pvalues[,2],"BH")),t='h', col=cols[map[,"Chr"]])                                                                 # BH corrected for multiple testing
#  abline(h=-log10(0.05), col="orange", lwd=1,lty=2) ; abline(h=-log10(0.01), col="green", lwd=1,lty=2)                                  # Significance thresholds
}

# Look into our top marker
plot(phenotype ~ as.factor(unlist(genotypes[which(-log10(pvalues[,2]) > 5),])))