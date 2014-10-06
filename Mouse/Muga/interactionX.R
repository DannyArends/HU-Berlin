# interaction QTLs
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

setwd("E:/Mouse/DNA/MegaMuga/")                                                                               # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE, colClasses="character")      # Normal A, H, B genotypes

setwd("E:/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI")                                                # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),phenos]                   # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)]                                                 # The F2 individuals
F2 <- F2[-which(F2=="6661459")]                                                                               # This individual has no genotype data

genotypes   <- genotypes[,F2]

setwd("E:/Mouse/ClassicalPhenotypes/AIL")
qtls <- read.table("Analysis/qtls_fatDlean_gwas.txt", sep="\t", colClasses=c("character",rep("numeric",3)), header=TRUE)

topChr3 <- rownames(qtls[which.max(qtls[,"qtl70"]),])
markersOnX <- rownames(map[(map[,"Chr"]=="X"),])

ind           <- colnames(genotypes[topChr3,!is.na(genotypes[topChr3,])])
phenotype     <- phenotypes[ind, paste0("mri42d_fat")] /  phenotypes[ind, paste0("mri70d_lean")]
M1            <- as.factor(as.character(genotypes[topChr3,ind]))
littersize    <- as.numeric(phenotypes[ind, "WG2"])

hetrozygousness <- apply(genotypes[markersOnX,],1,function(x){length(which(as.character(x)=="H"))})
markersOnXNoH <- markersOnX[which(!markersOnX %in% names(which(hetrozygousness > 20)))]

interactions  <- NULL
for(marker in markersOnXNoH){
  M2     <- as.factor(as.character(genotypes[marker, ind]))
  tryCatch(res  <- anova(lm(phenotype ~  littersize + M1 + M1:M2 ))[[5]][3], error = function(e){ res <<- NA })
  interactions <- c(interactions, -log10(res))
}
names(interactions) <- markersOnXNoH                   # Name the markers
topChrX <- interactions[which.max(interactions)]               # LOD score of 6
plot(interactions,t='h')

colorz <- c("red", "darkgreen", "blue")
names(colorz) <- c("A","H","B")

correctedphenotype <- lm(phenotype ~  littersize)$residuals

plot(c(1,3.5), c(-0.3, 0.5), t='n', xaxt="n")
for(elemX in c("A","H","B")){
  elemLine <- NULL
  moffset <- match(elemX, c("A","H","B"))/20
  for(elem in c("A","H","B")){
    ix <- which(genotypes[topChr3,] == elem & genotypes[names(topChrX),] == elemX)
    mmean <- mean(correctedphenotype[ix],na.rm=TRUE)
    msd <- sd(correctedphenotype[ix],na.rm=TRUE)
    cat(length(ix), elem, elemX, mmean,"\n")
    points(match(elem,c("A","H","B"))+moffset, mmean, col=colorz[elem])
    points(match(elem,c("A","H","B"))+moffset, mmean+msd, pch="-", col=colorz[elem])
    points(match(elem,c("A","H","B"))+moffset, mmean-msd, pch="-", col=colorz[elem])
    points(rep(match(elem,c("A","H","B")),3)+moffset, c(mmean-msd,mmean,mmean+msd),type='l', col=colorz[elem])
    elemLine <-c(elemLine, mmean)
    
    ix <- which(genotypes[topChr3,] == elem)
    mmean <- mean(correctedphenotype[ix],na.rm=TRUE)
    msd <- sd(correctedphenotype[ix],na.rm=TRUE)
    points(match(elem,c("A","H","B"))+moffset, mmean, col="black",pch="X")
  }
  points(c(1,2,3)+moffset, elemLine, type='l', col=colorz[elemX])
}

# Get the region
names(which(interactions> 4))
# Flank: backupUNC200098152   X  69213997 35.2874
# Top:         UNC200019054   X  74112481 37.7436
# Flank: backupUNC200099374   X  74266257 37.9316

### Get the genes in the region, using the snpToGene files (see snpToGene.R)
sel <- which(EXONS[which(EXONS[,1] == "X"),4] > 69213000 &  EXONS[which(EXONS[,1] == "X"), 5] < 74266300)
geneids <- unique(gsub("gene_id ","",unlist(lapply(strsplit(as.character(EXONS[which(EXONS[,1] == "X"),][sel, 9]),"; "),"[",1))))
RPKM[which(RPKM[,1] %in% geneids),"mgi_description"]


iii <- which(apply(RPKM[which(RPKM[,1] %in% geneids),c(35:40)],1,function(x){return(sum(as.numeric(x))); }) > 0)
RPKM[which(RPKM[,1] %in% geneids),c(7,8,35:40)][iii,]