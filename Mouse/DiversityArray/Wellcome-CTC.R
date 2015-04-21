# Analysis of Wellcome-CTC data from the mouse diversity array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#

setwd("E:/Mouse/DNA/DiversityArray/")
wellcomeannot <- read.table("CTC/RAW/final-9-2-5.csv", sep=",", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)

# Create a fasta file with the probe sequences to blast against the reference genome
if(!file.exists("CTC/Analysis/probes.fasta")){
  cat("", file="CTC/Analysis/probes.fasta")
  for(x in 1:nrow(wellcomeannot)){
    if(as.character(wellcomeannot[x,"Sequence"]) != ""){
      cat(">OLD", rownames(wellcomeannot)[x], "\n", sep = "", file="CTC/Analysis/probes.fasta", append=TRUE)
      probesequence <- unlist(strsplit(wellcomeannot[x, "Sequence"],""))
      if(length(which(probesequence == "["))){
        fsequence <- paste0(paste0(probesequence[1:(which(probesequence == "[")-1)], collapse=""), "N", paste0(probesequence[(which(probesequence == "]")+1):length(probesequence)], collapse=""))
        cat(fsequence, "\n", sep = "", file="CTC/Analysis/probes.fasta", append=TRUE)
      }
    }
  }
}

# blastn -task blastn -query CTC/Analysis/probes.fasta -db E:/Mouse/DNA/DiversityArray/Analysis/Mus_musculus.GRCm38.74.dna.db -perc_identity 95 -outfmt 6 -evalue=0.1 -out CTC/Analysis/probelocations.txt

# Load the locations of the probes and filter them so that we only keep the unique matching probes
locations <- read.csv("CTC/Analysis/probelocations.txt", sep = "\t", header=FALSE)
colnames(locations) <- c("ProbeName", "Chr", "Ident", "Length", "U1", "U2", "U3", "Match", "Start", "Stop", "evalue", "Score")

### BLAST post-processing
locations   <- locations[-which(locations[,"Score"] < 60),]                                                                                     # Match is not good enough to be considered as duplicate ( evalue < 60 )
dupprobes   <- unique(locations[which(duplicated(locations[,"ProbeName"])),"ProbeName"])                                                        # Probes which have multiple matches
goodprobes  <- locations[ which(!locations[,"ProbeName"] %in% dupprobes),]                                                                      # No multi-mapping probes
goodprobes  <- goodprobes[, c("ProbeName", "Chr", "Start", "Stop")]
rownames(goodprobes) <- goodprobes[,"ProbeName"] ; goodprobes <- goodprobes[,-1]                                                                # Use first column 'ProbeName' as rownames
goodprobes <- cbind(goodprobes, Location = floor((goodprobes[,"Start"] + goodprobes[,"Stop"])/2))

### Get the reference and alternative allele from the probe data
refalt <- matrix("", nrow(goodprobes), 2, dimnames=list(rownames(goodprobes), c("Ref","Alt")))
for(x in rownames(goodprobes)){
  probesequence <- unlist(strsplit(wellcomeannot[gsub("OLD","",x), "Sequence"],""))
  refalt[x, "Ref"] <- probesequence[(which(probesequence == "[")+1)] ; refalt[x, "Alt"] <- probesequence[(which(probesequence == "]")-1)]
}
goodprobes <- cbind(goodprobes, refalt)

if(!file.exists("CTC/Analysis/alldata.txt")){
  snpdata1 <- read.table("CTC/RAW/Strains-0.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
  snpdata2 <- read.table("CTC/RAW/Strains-200.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
  snpdata3 <- read.table("CTC/RAW/Strains-400.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
  wellcomedata <- cbind(snpdata1, snpdata2, snpdata3)[, -c(203, 204, 405, 406)]
  write.table(alldata, file="CTC/Analysis/alldata.txt", sep="\t")
}else{
  wellcomedata <- read.table(file="CTC/Analysis/alldata.txt", sep="\t", header = TRUE, row.names=1, colClasses="character", check.names=FALSE)
}
strains <- colnames(wellcomedata)[-c(1:2)]
cat("Wellcome data loaded:", length(strains), "strains,", nrow(wellcomedata), "SNPs\n")                                                         # Wellcome data loaded: 479 strains, 13368 SNPs
wellcomedata <- wellcomedata[rownames(wellcomedata) %in% gsub("OLD","",rownames(goodprobes)),]
cat("Wellcome data loaded:", length(strains), "strains,", nrow(wellcomedata), "SNPs after QC\n")                                                # Wellcome data loaded: 479 strains, 13368 SNPs
goodprobes <- goodprobes[gsub("OLD","",rownames(goodprobes)) %in% rownames(wellcomedata), ]

nwellcome <- matrix(NA, nrow(wellcomedata), ncol(wellcomedata) + 4, 
                    dimnames = list(rownames(wellcomedata), c("Chr", "Location", "CTC_A", "CTC_B", "Ref", "Alt", colnames(wellcomedata)[-c(1,2)])))

for(x in 1:nrow(wellcomedata)){
  name <- paste0("OLD", rownames(nwellcome)[x])

  nwellcome[x,"Chr"] <- as.character(goodprobes[name,"Chr"])
  nwellcome[x,"Location"] <- as.character(goodprobes[name,"Location"])
  gtype <- unique(unlist(wellcomedata[x, 3:ncol(wellcomedata)]))
  gtypes <- as.character(wellcomedata[x, 3:ncol(wellcomedata)])
  homo <- gtype[which(gtype != "H")]
  nwellcome[x, "CTC_A"] <- homo[1] ; nwellcome[x, "CTC_B"] <- homo[2]
  nwellcome[x, "Ref"] <- as.character(goodprobes[name,"Ref"]); nwellcome[x, "Alt"]  <- as.character(goodprobes[name, "Alt"])
  nwellcome[x, which(gtypes == homo[1]) + 6] <- 0
  nwellcome[x, which(gtypes == "H") + 6] <- 1
  nwellcome[x, which(gtypes == homo[2]) + 6] <- 2
  if((x %% 500) == 0) cat(x,"\n")
}
nwellcome <- data.frame(nwellcome, check.names = FALSE, stringsAsFactors = FALSE)
nwellcome[, 7:ncol(nwellcome)] <- apply(nwellcome[, 7:ncol(nwellcome)], 2, as.numeric)

### Create a dendrogram of all lines (100 random SNPs)
clusters   <- hclust(dist(t(nwellcome[sample(nrow(nwellcome), 100), strains])))                                                       # 100 random SNPs
dendrogram <- as.dendrogram(clusters)                                                                                                           # Dendrogram of all ALL lines
par(cex = 0.2, font=3)
plot(dendrogram)

### Create a file that can be used to get the rsIDs using the VEP webtool
if(!file.exists("CTC/Analysis/toAnnotate.txt")){
  cat("", file="CTC/Analysis/toAnnotate.txt")
  for(x in 1:nrow(goodprobes)){
    cat(as.character(goodprobes[x,"Chr"]), as.character(goodprobes[x,"Location"]), as.character(goodprobes[x,"Location"]), 
        paste0(goodprobes[x,"Ref"],"/",goodprobes[x,"Alt"]),1, "\n",file="CTC/Analysis/toAnnotate.txt", append=TRUE)
  }
}

### Use the VEP to annotate our file: http://www.ensembl.org/Mus_musculus/Tools/VEP
VEPout <- read.table("CTC/Analysis/VEPoutput.txt",sep="\t", skip=1,header = FALSE)

if(!file.exists("CTC/Analysis/alldata_annotated.txt")){
  wellcomedata <- cbind(nwellcome[,1:6], rsID = NA, nwellcome[,strains])
  rownames(wellcomedata) <- paste0("OLD", rownames(wellcomedata))

  for(x in 1:nrow(wellcomedata)){
    name <- paste(wellcomedata[x, "Chr"], wellcomedata[x, "Location"], paste0(wellcomedata[x, "Ref"],"/", wellcomedata[x, "Alt"]), sep="_")
    if(name %in% as.character(VEPout[,1])){
      #cat(name, as.character(VEPout[which(as.character(VEPout[,1]) == name)[1], 20]),"\n")
      wellcomedata[x,"rsID"] <- as.character(VEPout[which(as.character(VEPout[,1]) == name),20][1])
    }else{
      cat(name,"\n")
    }
  }
  write.table(wellcomedata, file="CTC/Analysis/alldata_annotated.txt", sep = "\t", quote = FALSE)
}else{
  wellcomedata <- read.table("CTC/Analysis/alldata_annotated.txt", sep="\t", row.names = 1, header=TRUE, check.names=FALSE)
}

### Calculate the genotype - genotype correlation matrix
if(!file.exists("CTC/Analysis/genoCorrelation.txt")){                                                                                           # Genotype correlations
  corM <- matrix(0, nrow(nwellcome), nrow(nwellcome), dimnames = list(rownames(wellcomedata), rownames(wellcomedata)))
  for(x in 1:nrow(nwellcome)){
    corV <- as.numeric(cor(as.numeric(nwellcome[x, 5:ncol(nwellcome)]), t(nwellcome[, 5:ncol(nwellcome)]), use="pair"))
    corM[x, ] <- corV
    if((x %% 50) == 0) cat(x, "\n")
  }
  write.table(round(corM, d = 2), "CTC/Analysis/genoCorrelation.txt", sep="\t", quote=FALSE)
}


### Load in our own data
ourdata <- read.table(file="Analysis/measurementsALL_annotated.txt", sep="\t", header=TRUE)

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                                                                         # Load the annotation
annInData  <- match(colnames(ourdata), rownames(annotation))                                                                                    # Align JAX B6 with the JAX BMMI data
annotation <- annotation[annInData,]
cat("Our data loaded:", ncol(ourdata), "strains,", nrow(ourdata), "SNPs\n")                                                                     # Our data loaded: 86 strains, 543047 SNPs

wellcomeoverlap <- wellcomedata[which(as.character(wellcomedata[,"rsID"]) %in% as.character(ourdata[,"dbSNP_ID"])),]
wellcomeoverlap <- wellcomeoverlap[!is.na(wellcomeoverlap[,"rsID"]), ]
wellcomeoverlap <- wellcomeoverlap[!duplicated(wellcomeoverlap[,"rsID"]), ]
rownames(wellcomeoverlap) <- as.character(wellcomeoverlap[,"rsID"])

ourdataoverlap  <- ourdata[which(as.character(ourdata[,"dbSNP_ID"]) %in% as.character(wellcomedata[,"rsID"])),]
ourdataoverlap  <- ourdataoverlap[!is.na(ourdataoverlap[,"dbSNP_ID"]), ]
ourdataoverlap  <- ourdataoverlap[!duplicated(ourdataoverlap[,"dbSNP_ID"]), ]
rownames(ourdataoverlap) <- as.character(ourdataoverlap[,"dbSNP_ID"])

cat("Our data loaded:", ncol(ourdataoverlap), "strains,", nrow(ourdataoverlap), "SNPs\n")                                                       # Our data loaded: 86 strains, 3335 SNPs
cat("Wellcome loaded:", ncol(wellcomeoverlap), "strains,", nrow(wellcomeoverlap), "SNPs\n")                                                     # Wellcome loaded: 86 strains, 3335 SNPs


### Create a plot
chromosomes  <- as.character(c(1:19, "X", "Y", "M"))
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)

plot(c(0, max(chrInfo[,"Length"])), c(1,nrow(chrInfo)), t='n', yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(wellcomeoverlap, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  points(x=xloc, y=(yloc-0.2), pch='|', col='black',cex=1)
})


