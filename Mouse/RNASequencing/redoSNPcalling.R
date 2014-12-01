fileBases <- c("5068_GAGTGG_L004_", "5069_AGTCAA_L004_",                                # BFMI samples
               "4868_GCCAAT_L001_", "5067_ATCACG_L004_",                                # B6N samples
               "5070_CGATGT_L005_", "5071_CCGTCC_L005_", "5072_TAGCTT_L005_",           # maternal B6N samples
               "5073_TTAGGC_L006_", "5074_GATCAG_L006_", "5075_ATGTCA_L006_")           # maternal BFMI samples

gatk          <- "/opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"
referenceDir  <- "genomes"
referenceName <- "Mus_musculus.GRCm38.74.dna"
reference     <- paste0(referenceDir, "/", referenceName, ".fasta")
knownsnps     <- paste0(referenceDir, "/mgp.v3.snps.rsIDdbSNPv137.vcf")                                # Reference SNPs, download from: ftp://ftp-mouse.sanger.ac.uk/
mpileup       <- "~/Github/samtools/samtools mpileup"

outputSIRBAMs <- NULL
for(fileBase in fileBases){
  outputSIRBAMs  <- c(outputSIRBAMs, paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam"))
}

cat(paste0("nohup ", mpileup, " -g -f ",reference," ", paste(outputSIRBAMs, collapse=" "), " -o population.bcf &"), "\n")
cat(paste0("nohup bcftools call -c population.bcf | ~/Github/bcftools/vcfutils.pl varFilter -d 10 - > population.vcf &"), "\n")


createNames <- function(x){ paste0(x[,1],":", x[,2],"_", x[,5]) }
getGenotypes <- function(x){ unlist(lapply(strsplit(x, ":"),"[",1)) }
getProbabilities <- function(x){ unlist(lapply(strsplit(x, ":"),"[",2)) }
getMaxProb <- function(x){ 
  p <- getProbabilities(x)
  unlist(lapply(strsplit(p, ","),function(i){max(as.numeric(i))}))
}

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/ReAnalysisSNPs")
vcfdata <- read.table("population.vcf", colClasses="character")
colnames(vcfdata) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","5068","5069","4868","5067","5070","5071","5072","5073","5074","5075")
rownames(vcfdata) <- createNames(vcfdata)
parents <- c("5068","5069","4868","5067")
matB6N  <- c("5070","5071","5072")
matBFMI <- c("5073","5074","5075")
samples <- c("5068","5069","4868","5067","5070","5071","5072","5073","5074","5075")

genos <- matrix(NA, nrow(vcfdata), length(samples))
probs <- matrix(NA, nrow(vcfdata), length(samples))
colnames(genos) <- samples ; rownames(genos) <- createNames(vcfdata)
colnames(probs) <- samples ; rownames(probs) <- createNames(vcfdata)

for(s in samples){ 
  genos[,s] <- getGenotypes(vcfdata[,s])
  probs[,s] <- getMaxProb(vcfdata[,s])
}

probInParents <- which(apply(probs[,parents], 1,function(x){ sum(x > 50) == 4 }))
genos <- genos[probInParents, ] ; probs <- probs[probInParents, ]

difGenoInParents <- which(apply(genos[,parents], 1,function(x){ 
  if(any(x == "0/1")) return(FALSE)
  if(x[1] == x[2] && x[3] == x[4] && x[1] != x[3]) return(TRUE)
  return(FALSE)
}))

genos <- genos[difGenoInParents, ] ; probs <- probs[difGenoInParents, ]

genos[which(probs < 50)] <- NA

showsASE <- apply(genos, 1, function(x){
  naB6N <- is.na(x[matB6N])
  naBFMI <- is.na(x[matBFMI])
  if(sum(naB6N) < 2){
    if(any(x[matB6N][!naB6N] == "0/1")) return(FALSE)
    return(TRUE)
  }
  if(sum(naBFMI) < 2){
    if(any(x[matBFMI][!naBFMI] == "0/1")) return(FALSE)
    return(TRUE)
  }
  return(FALSE)
})

genos <- genos[unlist(showsASE), ]

ASE <- t(apply(genos,1,function(x){
  naB6N <- is.na(x[matB6N])
  naBFMI <- is.na(x[matBFMI])
  matB6Ngeno <- "-"
  matBFMIgeno <- "-"
  if(sum(naB6N) < 2){
    tbl <- table(x[matB6N][!naB6N])
    if(length(tbl) > 1){
      matB6Ngeno <- "MIX"
    }else{
      if(x[1] == names(tbl)) matB6Ngeno <- "B6N"
      if(x[3] == names(tbl)) matB6Ngeno <- "BFMI"
    }
  }
  if(sum(naBFMI) < 2){
    tbl <- table(x[matBFMI][!naBFMI])
    if(length(tbl) > 1){
      matBFMIgeno <- "MIX"
    }else{
      if(x[1] == names(tbl)) matBFMIgeno <- "B6N"
      if(x[3] == names(tbl)) matBFMIgeno <- "BFMI"
    }
  }
  return(c(matB6Ngeno, matBFMIgeno))
}))
colnames(ASE) <- c("matB6N","matBFMI")

genos <- cbind(vcfdata[rownames(genos),c("CHROM","POS","REF","ALT")], genos, ASE)

hasASE <- which(apply(genos[,c("matB6N", "matBFMI")],1,function(x){
  any(x == "BFMI" | x == "B6N")
}))

genos <- genos[hasASE,]

setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
RPKM        <- read.csv("Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")   # RPKM values from RNA-Seq

GTF <- read.table("GTF/Mus_musculus.GRCm38.76.gtf", sep="\t")                                                       # Gene models
EXONS <- GTF[which(GTF[,3]=="exon"),]

datastr <- strsplit(as.character(EXONS[,9]), "; ")
lengths <- unlist(lapply(datastr, length))

addInfo <- matrix(NA, length(datastr), 2)
for(x in 1:length(datastr)){ addInfo[x, ] <- c(strsplit(datastr[[x]][1]," ")[[1]][2], strsplit(datastr[[x]][lengths[x]]," ")[[1]][2]); }

EXONS <- cbind(addInfo, EXONS)
uniqueGenes <- as.character(unique(EXONS[,1]))

geneExonsCoupling <- vector("list", length(uniqueGenes))                                                        # Create a list which enumerate the exons per gene
x <- 1
for(gene in uniqueGenes){
  geneExons <- which(EXONS[,1] == gene)
  possibleExons <- EXONS[geneExons[!duplicated(EXONS[geneExons,2])], -11]
  geneExonsCoupling[[x]] <- possibleExons
  cat(x, "Gene", gene, "has", nrow(geneExonsCoupling[[x]]), "/", length(geneExons),"Exons\n")
  x <- x + 1
}

geneSNPCoupling <- vector("list", length(uniqueGenes))
x <- 1
for(gene in geneExonsCoupling){                                                                               # Summarize all SNPs from a direction per gene
  onChr <- as.character(genos[,"CHROM"]) == as.character(unique(gene[,3]))
  for(exon in 1:nrow(gene)){
    exonStart <- as.numeric(gene[exon,6])
    exonEnd <- as.numeric(gene[exon,7])
    inEXON <- which(genos[onChr,"POS"] >= exonStart & genos[onChr,"POS"] <= exonEnd)
    if(length(inEXON) > 0){
      geneSNPCoupling[[x]] <- rbind(geneSNPCoupling[[x]], cbind(snpID = rownames(genos[onChr,][inEXON,]), genos[onChr,][inEXON,], Exon=exon))
      geneSNPCoupling[[x]] <- geneSNPCoupling[[x]][!duplicated(geneSNPCoupling[[x]][,"snpID"]),]
    }
    #names(geneSNPCoupling)[x] <- as.character(unique(gene[,1]))
  }
  if(!is.null(geneSNPCoupling[[x]])) cat(x, "found", nrow(geneSNPCoupling[[x]]), "SNPs in gene\n")
  x <- x + 1
}

names(geneSNPCoupling) <- uniqueGenes

geneSNPshort <- NULL
for(x in 1:length(geneSNPCoupling)){
  if(!is.null(geneSNPCoupling[[x]])) geneSNPshort <- c(geneSNPshort, geneSNPCoupling[x])
}

consistentASE <- which(unlist(lapply(geneSNPshort, function(x){
  if(length(table(as.character(x[,"matB6N"]))) == 1 || length(table(as.character(x[,"matBFMI"]))) == 1) return(TRUE)
  return(FALSE)
})))

geneSNPshort <- geneSNPshort[consistentASE]

matrixform <- NULL
for(x in 1:length(geneSNPshort)){
  additionalgeneInfo <- RPKM[which(RPKM[,"ensembl_gene_id"] == names(geneSNPshort)[x]), -c(1)]
  matrixform <- rbind(matrixform, cbind(ensembl_gene_id = names(geneSNPshort)[x], cbind(geneSNPshort[[x]]), additionalgeneInfo))
}

matrixform[,samples] <- apply(matrixform[,samples],2,function(x){gsub("/", "|", x)})
write.table(matrixform,"RPKM+ASE.txt", sep="\t", quote=FALSE,row.names=FALSE)

# ASE plots
chromosomes  <- as.character(c(1:19, "X", "Y", "MT"))

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo  <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
markers  <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
mlength  <- max(chrInfo[,"Length"])

plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="Allele Specific Expression", sub="maternal B6N (left, B6N♀ x BFMI♂) versus maternal BFMI (right, BFMI♀ x B6N♂)", yaxt="n", xlab="", ylab="", xaxt="n")
abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

aa <- apply(matrixform, 1,function(x){
  xloc <- match(as.character(x["CHROM"]), chromosomes); yloc <- as.numeric(x["POS"])
    col <- "white"
    if(x["matB6N"]=="B6N")  col <- "gray60"
    if(x["matB6N"]=="BFMI") col <- "gold1"
    if(col != "white") points(x=xloc-0.2, y=yloc, pch="-", col=col, cex=2.0)
    col <- "white"
    if(x["matBFMI"]=="B6N")  col <- "gray60"
    if(x["matBFMI"]=="BFMI") col <- "gold1"
    if(col != "white") points(x=xloc+0.2, y=yloc, pch="-", col=col, cex=2.0)
})

cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="black", lty=1, lwd=2)
  cnt <<- cnt + 1
})

axis(1,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1, cex.axis=1.5)
axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=1.2, las=1)
legend("topright", c("BFMI allele expressed", "B6N allele expressed"), fill=c("gold1","gray60"), cex=1.2)

library(biomaRt)                                                                                        # Biomart package
library(topGO)                                                                                          # topGO package
  
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
allgenes <- RPKM[, "ensembl_gene_id"]

if(!file.exists("GeneOntology/GOannotation.txt")){
  bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                              # Biomart for mouse genes
  biomartResults <- NULL
  for(x in seq(1, length(allgenes), 1000)){                                                             # Do 1000 per time, just to please biomaRt
    xend <- min((x + 1000),length(allgenes))                                                            # Don't walk passed the end of the array
    cat("Retrieving", x, "/", xend,"\n")
    res.biomart <- getBM(c("ensembl_gene_id","go_id"),                                                  # Use biomart to retrieve GO terms
                          filters="ensembl_gene_id", values=allgenes[x:xend], mart=bio.mart)
    biomartResults <- rbind(biomartResults, res.biomart)
  }
  write.table(biomartResults, file="GeneOntology/GOannotation.txt", sep="\t", row.names=FALSE)
}else{
  cat("Loading biomaRt gene ontology annotation from disk\n")
  biomartResults <- read.table("GeneOntology/GOannotation.txt", sep="\t", header=TRUE)
}

if(!file.exists("GeneOntology/geneid2go.map")){                                                                      # Create the ENSG to GO map only if it doesn't exists
  cat("", file="GeneOntology/geneid2go.map")
  for(ensid in unique(biomartResults[,"ensembl_gene_id"])){
    idxes <- which(biomartResults[,"ensembl_gene_id"] == ensid)
    goids <- biomartResults[idxes,"go_id"]
    emptygo <- which(goids=="")
    if(length(emptygo) > 0) goids <- goids[-emptygo]
    if(length(goids) > 0) cat(ensid,"\t", paste(goids, collapse=", "),"\n", file="GeneOntology/geneid2go.map", append=TRUE, sep="")
  }
}

# Do gene ontology
doGO <- function(allgenes, selected){
  genelist <- rep(0, length(allgenes))                                                                # Create a gene list
  names(genelist) <- allgenes                                                                         # Add the names
  genelist[selected] <- 1                                                                             # Set the switched genes to 1

  geneID2GO     <- readMappings(file = "GeneOntology/geneid2go.map")
  GOdata        <- new("topGOdata", ontology = "BP", allGenes = as.factor(genelist), annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultFisher  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

  #pdf("GeneOntologyTree.pdf")
    showSigOfNodes(GOdata, topGO::score(resultFisher), firstSigNodes = 5, useInfo = 'all')
  #dev.off()
  return(list(GOdata,resultFisher))
}

ASEgenesGO <- doGO(allgenes, unique(matrixform[,"ensembl_gene_id"]))
ASEgenes   <- GenTable(ASEgenesGO[[1]], classicFisher = ASEgenesGO[[2]], orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 20)

maternal <- which(matrixform[,"matB6N"] == "B6N" &  matrixform[,"matBFMI"] == "BFMI")
ASEmaternalGO <- doGO(allgenes, unique(matrixform[maternal, "ensembl_gene_id"]))
ASEmaternal   <- GenTable(ASEmaternalGO[[1]], classicFisher = ASEmaternalGO[[2]], orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 20)

paternal <- which(matrixform[,"matB6N"] == "BFMI" &  matrixform[,"matBFMI"] == "B6N")
ASEpaternalGO <- doGO(allgenes, unique(matrixform[paternal, "ensembl_gene_id"]))
ASEpaternal   <- GenTable(ASEpaternalGO[[1]], classicFisher = ASEpaternalGO[[2]], orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 20)

alwaysBFMI <- which(matrixform[,"matB6N"] == "BFMI" &  matrixform[,"matBFMI"] == "BFMI")
ASEbfmiGO <- doGO(allgenes, unique(matrixform[alwaysBFMI, "ensembl_gene_id"]))
ASEbfmi   <- GenTable(ASEbfmiGO[[1]], classicFisher = ASEbfmiGO[[2]], orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 20)

alwaysB6N <- which(matrixform[,"matB6N"] == "B6N" &  matrixform[,"matBFMI"] == "B6N")
ASEb6nGO <- doGO(allgenes, unique(matrixform[alwaysB6N, "ensembl_gene_id"]))
ASEb6n   <- GenTable(ASEb6nGO[[1]], classicFisher = ASEb6nGO[[2]], orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 20)

