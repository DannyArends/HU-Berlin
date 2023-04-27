setwd("/home/rqdt9/Data/BFMI")
genesmodelgff <- read.csv("Mus_musculus.GRCm38.89.chromosome.3.gff3", comment.char="#", sep="\t", colClasses="character", header=FALSE)
#genesmodelgff <- genesmodelgff[genesmodelgff[,"V3"] %in% c("gene", "mRNA", "miRNA","exon", "CDS", "five_prime_UTR", "three_prime_UTR"),]
genesmodelgff <- genesmodelgff[sort(as.numeric(genesmodelgff[,"V4"]), index.return=TRUE)$ix,]

tfbsgff <- read.csv("mus_musculus.GRCm38.motiffeatures.20161111.gff", header=FALSE, sep="\t")
tfbsgff <- tfbsgff[which(tfbsgff[,"V1"] == 3),]

reggff <- read.csv("mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20161111.gff", header=FALSE, sep="\t")
reggff <- reggff[which(reggff[,"V1"] == 3),]

getRegion <- function(genesmodelgff, pos.start, pos.stop){
  region <- genesmodelgff[which(as.numeric(genesmodelgff[,"V4"]) >= pos.start & as.numeric(genesmodelgff[,"V5"]) <=  pos.stop),]
  extraD <- unique(unlist(lapply(strsplit(region[,"V9"],";"), function(x){unlist(lapply(strsplit(x,"="),"[",1))})))
  extraColumns <- matrix(NA, nrow(region), length(extraD), dimnames=list(rownames(region), extraD))
  region <- cbind(region, extraColumns)

  # Unpack the extra column into it's own columns, so we can more easily filter on the data
  for(x in 1:nrow(region)){
    extras <- unlist(strsplit(region[x,"V9"],";"))
    esplit <- strsplit(extras, "=")
    keys <- unlist(lapply(esplit,"[",1))
    values <- unlist(lapply(esplit,"[",2))
    region[x, keys] <- values
    if(x %% 1000 == 0) cat(x, "/", nrow(region), "\n")
  }
  region <- region[,-which(colnames(region) == "V9")]
}
substrRight <- function(x, n) { substr(x, nchar(x)-n+1, nchar(x)) }

region.fm <- getRegion(genesmodelgff, 36480000, 36850000)

genesInRegion <- function(region, includeGM = FALSE, includeRIKEN = FALSE, verbose=FALSE){
  allgenes <- region[which(region[,"V3"] == "gene"), "Name"]
  if(!includeGM){
    gmGenes <- which(substr(allgenes,1,2) == "Gm" & nchar(allgenes) >= 6)
    if(verbose) cat("Removing", length(gmGenes), "Gm genes\n")
    if(length(gmGenes) >0) allgenes <- allgenes[-gmGenes]
  }
  if(!includeRIKEN){
    rikenGenes <- which(substrRight(allgenes,3) == "Rik")
    if(verbose) cat("Removing", length(rikenGenes), "riken genes\n")
    if(length(rikenGenes) >0) allgenes <- allgenes[-rikenGenes]
  }
  
  genedata <- NULL
  for(gene in allgenes) {
    genedata <- rbind(genedata, region[which(region[,"Name"] == gene),])
  }
  #emptycols <- apply(genedata,2,function(x){sum(is.na(x)) == length(x)})
  #genedata <- genedata[,-which(emptycols)]
  return(genedata)
}

genes.fm <- genesInRegion(region.fm)

getTranscripts <- function(region, gene.name = "Bbs7"){
  genes <- genesInRegion(region)
  idx <- which(genes[,"Name"] == gene.name)
  if(length(idx) != 1) stop(paste0("Gene '",gene.name,"' not found in region"))
  inG <- region[which(region[,"Parent"] == genes[idx,"ID"]),]
  inG <- inG[which(inG[, "V3"] == "mRNA"),]
  transcriptIDs <- gsub("transcript:","", inG[, "ID"]) #region[which(region[,"Parent"] == genes[idx,"ID"]),"ID"]
  transcripts <- unique(na.omit(transcriptIDs))
  return(transcripts)
}

bbs7.transcripts <- getTranscripts(region.fm, "Bbs7")
ccna2.transcripts <- getTranscripts(region.fm, "Ccna2")

getTranscriptData <- function(region, transcript, what = "exon") {
  transcript.data <- region[grep(transcript, region[,"Parent"]),]
  transcript.data <- transcript.data[which(transcript.data[,"V3"] == what),]
  return(transcript.data)
}

genes.start <- min(c(as.numeric(genes.fm[,"V4"]), as.numeric(genes.fm[,"V5"]))) #36598403 #min(c(as.numeric(genes.fm[,"V4"]), as.numeric(genes.fm[,"V5"])))
genes.stop <- max(c(as.numeric(genes.fm[,"V4"]), as.numeric(genes.fm[,"V5"]))) #36602354 #max(c(as.numeric(genes.fm[,"V4"]), as.numeric(genes.fm[,"V5"])))

tfbsgff <- tfbsgff[which(tfbsgff[,"V4"] >= genes.start & tfbsgff[,"V5"] <= genes.stop),]
reggff <- reggff[which(reggff[,"V4"] >= genes.start & reggff[,"V5"] <= genes.stop),]

### load in the vcf file
vcf <- read.table("output.bbs7_33mb_39mb_region.vcf",sep="\t",colClasses="character", header = TRUE)
samples <- colnames(vcf)[-c(1:9)]

vcf <- vcf[which(as.numeric(vcf[,"QUAL"]) > 100),]  # High quality SNPs
vcf <- vcf[-grep(",", vcf[,"ALT"]),]         # only 1 alternative allele
vcf <- vcf[, -c(6, 7, 8, 9)]

del.from <- 36599483
del.to <- 36602058
bfmisamples <- samples[grep("BFMI", samples)]

vcft <- vcf[vcf[,"POS"] >= del.from & vcf[,"POS"] < del.to, ]

vcf[,samples] <- apply(vcf[,samples], 2,function(x) {
  unlist(lapply(strsplit(as.character(x), ":"), "[", 1))
})

allEqual <- apply(vcf[,samples],1,function(x) { all(x == x[1]) })
#vcf <- vcf[-which(allEqual),]

equalinBFMIs <- apply(vcf[,bfmisamples],1,function(x) { x[1] == x[2] && x[1] == x[3] })
#vcf <- vcf[equalinBFMIs,]

vcf <- cbind(vcf, inS12 = NA)
vcf[,"inS12"] <- vcf[,"BFMI860.S12"] == "1/1"


vcf <- cbind(vcf, inS1 = NA)
vcf[,"inS1"] <- vcf[,"BFMI860.S12"] == vcf[,"BFMI861.S1"]

vcf <- cbind(vcf, inS2 = NA)
vcf[,"inS2"] <- vcf[,"BFMI860.S12"] == vcf[,"BFMI861.S2"]


vcf <- cbind(vcf, inAKR_J = NA)
vcf <- cbind(vcf, AKR_JisBFMI = FALSE)
vcf[,"inAKR_J"] <- vcf[,"AKR_J"] == "1/1"#vcf[,bfmisamples[1]]
vcf[which(vcf[,"inAKR_J"]),"AKR_JisBFMI"] <- vcf[which(vcf[,"inAKR_J"]),"AKR_J"] == vcf[which(vcf[,"inAKR_J"]),bfmisamples[1]]

vcf <- cbind(vcf, inDBA_2J = NA)
vcf <- cbind(vcf, DBA_2JisBFMI = FALSE)
vcf[,"inDBA_2J"] <- vcf[,"DBA_2J"] == "1/1"#vcf[,bfmisamples[1]]
vcf[which(vcf[,"inDBA_2J"]),"DBA_2JisBFMI"] <- vcf[which(vcf[,"inDBA_2J"]),"DBA_2J"] == vcf[which(vcf[,"inDBA_2J"]),bfmisamples[1]]

vcf <- cbind(vcf, inNZO = NA)
vcf <- cbind(vcf, NZOisBFMI = FALSE)
vcf[,"inNZO"] <- vcf[,"NZO"] == "1/1"#vcf[,bfmisamples[1]]
vcf[which(vcf[,"inNZO"]),"NZOisBFMI"] <- vcf[which(vcf[,"inNZO"]),"NZO"] == vcf[which(vcf[,"inNZO"]),bfmisamples[1]]

vcf <- cbind(vcf, inSJL = NA)
vcf <- cbind(vcf, SJLisBFMI = FALSE)
vcf[,"inSJL"] <- vcf[,"SJL"] == "1/1"#vcf[,bfmisamples[1]]
vcf[which(vcf[,"inSJL"]),"SJLisBFMI"] <- vcf[which(vcf[,"inSJL"]),"SJL"] == vcf[which(vcf[,"inSJL"]),bfmisamples[1]]


vcf <- vcf[as.numeric(vcf[,"POS"]) > genes.start & as.numeric(vcf[,"POS"]) < genes.stop,]

against <- "AKR_J"

vcf[vcf[,2] > 36613477 & vcf[,2] < 36615477 ,]


vcf[vcf[,2] > 36589000 & vcf[,2] < 36590000 ,]

colorz <- RColorBrewer::brewer.pal(11, "Set3")
op <- par(cex=1)
op <- par(mar=c(5,8,5,5))
#intext <- paste0("BFMI complementation against ", against)
#op <- par(mar=c(4, 0.2, 0.2, 0.2))
#plot(c(36574142, 36618477), y = c(-0.2, 0.3), t ='n',yaxt = 'n', xaxt='n', ylab = "", xlab="Position (Mb)")
plot(c(36610000, 36618477), y = c(-0.2, 0.3), t ='n',yaxt = 'n', xaxt='n', ylab = "", xlab="Position (Mb)")
#plot(c(36598000, 36607500), y = c(-0.25, 0.3), t ='n',yaxt = 'n', xaxt='n', ylab = "", xlab="Position (Mb)")
#plot(c(36612477, 36615477), y = c(-0.2, 0.7), t ='n',yaxt = 'n', xaxt='n', ylab = "", xlab="Position (Mb)")
#plot(c(36618477, 36618477), y = c(-0.2, 0.7), t ='n',yaxt = 'n', xaxt='n', ylab = "", xlab="Position (Mb)")

#text(x = (genes.stop+genes.start)/2, y = 1.2, intext, cex=1.3)
#abline(h=0.5)
for(x in 1:nrow(genes.fm)){
  y <- (as.numeric(genes.fm[x, "V7"] == "+") * 1.2) - 0.1
  text(x = (as.numeric(genes.fm[x, "V4"]) + as.numeric(genes.fm[x, "V5"])) / 2 , y = 0.5 + 1.1 * (y - 0.5), genes.fm[x, "Name"], col=colorz[(x + 2)])
  transcripts <- getTranscripts(region.fm, genes.fm[x, "Name"])
  for(transcript in transcripts[1]){
    segments(as.numeric(genes.fm[x, "V4"]), y, as.numeric(genes.fm[x, "V5"]), lwd = 1)
    exons <- getTranscriptData(region.fm, transcript)
    for(e in 1:nrow(exons)){
      segments(as.numeric(exons[e, "V4"]), y, as.numeric(exons[e, "V4"]), lwd = 10, col=colorz[(x + 2)])
    }
    
    CDSs <- getTranscriptData(region.fm, transcript, "CDS")
    cat(transcript, nrow(CDSs), "\n")
    if(nrow(CDSs) > 0){
      for(cds in 1:nrow(CDSs)){
        segments(as.numeric(CDSs[cds, "V4"]), y, as.numeric(CDSs[cds, "V4"]), lwd = 4, col="gray")
      }
    }
    if(y > 0.5) y <- y - 0.1
    if(y < 0.5) y <- y + 0.1
    #Sys.sleep(1)
  }
}

vcfBFMI <- vcf[which(vcf[,bfmisamples[1]] != "./."), ]
#points(x = as.numeric(vcfBFMI[, "POS"]), y = rep(0.40, nrow(vcfBFMI)), pch = "|", cex=0.8, col="black")
points(x = as.numeric(vcf[vcf[,"inS12"],"POS"]), y = rep(0.0, length(which(vcf[,"inS12"]))), pch = "|", cex=0.8, col="black")
points(x = c(36599727, 36601264), y = c(0,0), t='l', col = "red")
points(x = c(36599890, 36601264), y = c(0.1,0.1), t='l', col = "red")
points(x = c(36600576, 36601264), y = c(0.05,0.05), t='l', col = "red")


vcfIN <- vcf[which(vcf[,paste0("in",against)]),]
vcfIN <- vcfIN[which(vcfIN[,against] != "./."), ]
#points(x = as.numeric(vcfIN[, "POS"]), y = rep(0.55, nrow(vcfIN)), pch = "|", cex=0.8, col="green")

vcfOUT <- vcf[which(!vcf[,paste0("in",against)]),]
vcfOUT <- vcfOUT[which(vcfOUT[,against] != "./."), ]
#points(x = as.numeric(vcfOUT[, "POS"]), y = rep(0.6, nrow(vcfOUT)), pch = "|", cex=0.8, col="orange")

points(x = as.numeric(vcf[vcf[,"inAKR_J"],"POS"]), y = rep(0.1, length(which(vcf[,"inAKR_J"]))), pch = "|", col=c("darkgreen","gray")[as.numeric(vcf[vcf[,"inAKR_J"],"AKR_JisBFMI"]) + 1], cex=0.7)
points(x = as.numeric(vcf[vcf[,"inDBA_2J"],"POS"]), y = rep(0.15, length(which(vcf[,"inDBA_2J"]))), pch = "|", col=c("blue","gray")[as.numeric(vcf[vcf[,"inDBA_2J"],"DBA_2JisBFMI"]) + 1], cex=0.7)
points(x = as.numeric(vcf[vcf[,"inNZO"],"POS"]), y = rep(0.05, length(which(vcf[,"inNZO"]))), pch = "|", col=c("orange","gray")[as.numeric(vcf[vcf[,"inNZO"],"NZOisBFMI"]) + 1], cex=0.7)
points(x = as.numeric(vcf[vcf[,"inSJL"],"POS"]), y = rep(0.2, length(which(vcf[,"inSJL"]))), pch = "|", col=c("blueviolet","gray")[as.numeric(vcf[vcf[,"inSJL"],"SJLisBFMI"]) + 1], cex=0.7)

legend("bottomright", rev(c( "BFMI", "NZO", "AKR/J", "DBA/2J", "SJL")), fill=rev(c("black", "orange", "darkgreen", "blue", "blueviolet")))

tfbst <- as.factor(unlist(lapply(strsplit(as.character(tfbsgff[,"V9"]), "="), "[", 3)))
regst <- as.factor(unlist(lapply(strsplit(as.character(reggff[,"V9"]), "="), "[", 6)))

#points(x = (as.numeric(tfbsgff[,"V4"]) + as.numeric(tfbsgff[,"V5"])) / 2, y = c(0.65, 0.35)[as.numeric(tfbsgff[,"V7"] == "+") + 1], pch = c("▼", "▲")[as.numeric(tfbsgff[,"V7"] == "+") + 1], col="black")
#points(x = (as.numeric(tfbsgff[,"V4"]) + as.numeric(tfbsgff[,"V5"])) / 2, y = c(0.65, 0.35)[as.numeric(tfbsgff[,"V7"] == "+") + 1], pch = c("▼", "▲")[as.numeric(tfbsgff[,"V7"] == "+") + 1], col=colorz[as.numeric(tfbst)], cex=0.8)
#legend("topright", levels(tfbst), pch="▼", col=colorz[1:nlevels(tfbst)], cex=0.8, pt.cex = 0.8)

#points(x = (as.numeric(reggff[,"V4"]) + as.numeric(reggff[,"V5"])) / 2, y = rep(0.7, nrow(reggff)), pch = 19, col=colorz[3+as.numeric(regst)])
#legend("bottomleft", levels(regst), pch=19, col=colorz[3+(1:nlevels(regst))], cex=0.8, pt.cex = 0.8)
#axis(2, at=c(0.5,0.55,0.6), c("BFMI",paste0(against, " == BFMI"),paste0(against, " != BFMI")), las=2, cex.axis=1)
xlocs <- seq(36550000, 36690000, 5000)
axis(1, at=xlocs, round(xlocs / 1000000, 2), las=1, cex.axis=1)

abline(v=36701498)
### Plot the different proteins

bbs7.transcripts <- getTranscripts(region.fm, "Bbs7")
ccna2.transcripts <- getTranscripts(region.fm, "Ccna2")
exocs9.transcripts <- getTranscripts(region.fm, "Exosc9")
trpc3.transcripts <- getTranscripts(region.fm, "Trpc3")

alltranscripts <- c(bbs7.transcripts, ccna2.transcripts, exocs9.transcripts, trpc3.transcripts)

for(transcript in alltranscripts){
  exons <- getTranscriptData(region.fm, transcript)
  cds <- getTranscriptData(region.fm, transcript, "CDS")
  ordering <- 1:nrow(cds)
  cds.start <- "V4"
  cds.end <- "V5"
  if(cds[1,"V7"] == "-"){
    ordering <- nrow(cds):1
    cds.start <- "V5"
    cds.end <- "V4"
  }
  
  bp <- 0
  for(x in ordering){
    # If the gene is on the negative strand the position is counted from the length of the CDS
    if(cds[1,"V7"] == "+"){
      bp <- bp + (1 + (as.numeric(cds[x,"V5"]) - as.numeric(cds[x,"V4"])))
    }else{
      bp <- bp + (1 + (as.numeric(cds[x,"V5"]) - as.numeric(cds[x,"V4"])))
    }
  }
  cat(transcript, cds[1,"V7"], length(ordering), bp, bp /3, "\n")
}


#primerlocations <- read.table("locations.txt")
#for(i in 1:nrow(primerlocations)){
#  points(points(primerlocations[i, "V3"], 0.7))
#  text(primerlocations[i, "V3"], 0.8, as.character(primerlocations[i, "V1"]))
#}
