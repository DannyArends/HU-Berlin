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

# blastn -task blastn -query CTC/Analysis/probes.fasta -db E:/Mouse/DNA/DiversityArray/Analysis/Mus_musculus.GRCm38.74.dna.db -perc_identity 100 -outfmt 6 -evalue=0.1 -out CTC/Analysis/probelocations.txt

# Load the locations of the probes and filter them so that we only keep the unique matching probes
locations <- read.csv("CTC/Analysis/probelocations.txt", sep = "\t", header=FALSE)
colnames(locations) <- c("ProbeName", "Chr", "Ident", "Length", "U1", "U2", "U3", "Match", "Start", "Stop", "evalue", "Score")

locations <- locations[-which(locations[,"Score"] < 80),]                                                                           # Match is not good enough to be considered as duplicate ( evalue < 80 )

dupprobes <- unique(locations[which(duplicated(locations[,"ProbeName"])),"ProbeName"])                                              # Probes which have multiple matches
bestlocs <- locations[ which(!locations[,"ProbeName"] %in% dupprobes),]                                                             # No multi-mapping probes

if(!file.exists("CTC/Analysis/alldata.txt")){
  snpdata1 <- read.table("CTC/RAW/Strains-0.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
  snpdata2 <- read.table("CTC/RAW/Strains-200.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
  snpdata3 <- read.table("CTC/RAW/Strains-400.09052005.txt", sep=" ", header=TRUE, colClasses="character", row.names=1, check.names=FALSE)
  wellcomedata <- cbind(snpdata1, snpdata2, snpdata3)[, -c(203, 204, 405, 406)]
  write.table(alldata, file="CTC/Analysis/alldata.txt", sep="\t")
}else{
  wellcomedata <- read.table(file="CTC/Analysis/alldata.txt", sep="\t", header = TRUE, row.names=1, colClasses="character", check.names=FALSE)
}
cat("Wellcome data loaded:", ncol(wellcomedata), "strains,", nrow(wellcomedata), "SNPs\n")                                          # Wellcome data loaded: 481 strains, 13368 SNPs

nwellcome <- matrix(NA, nrow(wellcomedata), ncol(wellcomedata) + 2, 
                    dimnames = list(rownames(wellcomedata), c("Chr", "Location", "CTC_A", "CTC_B", colnames(wellcomedata)[-c(1,2)])))

for(x in 1:nrow(wellcomedata)){
  nwellcome[x,"Chr"] <- as.character(wellcomedata[x,"CHR"])
  nwellcome[x,"Location"] <- as.character(wellcomedata[x,"POSITION"])
  gtype <- unique(unlist(wellcomedata[x, 3:ncol(wellcomedata)]))
  homo <- gtype[which(gtype != "H")]
  nwellcome[x, "CTC_A"] <- homo[1] ; nwellcome[x, "CTC_B"] <- homo[2]
  nwellcome[x, which(as.character(wellcomedata[x, 3:ncol(wellcomedata)]) == nwellcome[x, "CTC_A"]) + 4] <- 0
  nwellcome[x, which(as.character(wellcomedata[x, 3:ncol(wellcomedata)]) == "H") + 4] <- 1
  nwellcome[x, which(as.character(wellcomedata[x, 3:ncol(wellcomedata)]) == nwellcome[x, "CTC_B"]) + 4] <- 2
}
nwellcome <- data.frame(nwellcome, check.names = FALSE)


SNPregions <- apply(wellcomedata[1:10, c("CHR","POSITION","POSITION")], 1, paste0,collapse=":")

library(biomaRt)                                                                                                                    # Biomart
snp.db <- useMart("snp", dataset="mmusculus_snp")                                                                                   # For mouse SNPs


ourdata <- read.table(file="Analysis/measurementsALL_annotated.txt", sep="\t", header=TRUE,nrow=100)

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                          # Load the annotation
annInData  <- match(colnames(ourdata), rownames(annotation))                                     # Align JAX B6 with the JAX BMMI data
annotation <- annotation[annInData,]
cat("Our data loaded:", ncol(ourdata), "strains,", nrow(ourdata), "SNPs\n")                      # Our data loaded: 86 strains, 543047 SNPs

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)

plot(c(0, max(chrInfo[,"Length"])), c(1,nrow(chrInfo)), t='n', yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(wellcomedata, 1,function(x){
  yloc <- match(x["CHR"], chromosomes); xloc <- x["POSITION"]
  points(x=xloc, y=(yloc-0.2), pch='|', col='black',cex=0.4)
})


