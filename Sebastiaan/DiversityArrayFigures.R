# Figures created from the diversity array data for Sebastiaan
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Mar, 2015
# first written Mar, 2015
library(ape)

setwd("E:/Mouse/DNA/DiversityArray/")

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                           # Load the annotation

selected_lines <- c("X86022339", "X8609025", "X8613087", "X8613106", "X8523239", "X8563243", "C57BL.6CR_m_ep1")

rawdata <- read.table(file="Analysis/measurementsALL_annotated.txt", sep="\t", header=TRUE)
map      <- rawdata[,1:9]
calldata <- rawdata[,selected_lines]

rownames(map) <- map[,1]
rownames(calldata) <- map[,1]

# Use the Line names, not the individual IDs (Destructive do only once)
colnames(calldata) <- annotation[colnames(calldata), "Line"]

## TODO save so we do not have to reload the shit
 setwd("D:/Sebastiaan/")

clusters   <- hclust(dist(t(calldata), method = "manhattan"))
dendrogram <- as.dendrogram(clusters)                                                             # Dendrogram of all ALL lines
png("Dendrogram_SelectedLines.png", 1024, 768)
plot(dendrogram)
dev.off()
phylogram  <- as.phylo(clusters)                                                                  # Phylogram of all ALL lines
png("Cladogram_SelectedLines.png", 1024, 768)
plot(phylogram, type = "cladogram")
dev.off()

ordering <- rev(c("BFMI852", "BFMI856", "BFMI860-12", "BFMI860-S2", "BFMI861-S2", "BFMI861-S1"))
calldata <- calldata[,ordering]

hetro <- apply(calldata, 1, function(x){
  sum(x == 1, na.rm=TRUE)
})
map      <- map[-which(hetro >= 3),]
calldata <- calldata[-which(hetro >= 3),]

nas <- apply(calldata, 1, function(x){
  sum(is.na(x))
})
map      <- map[-which(nas >= 3),]
calldata <- calldata[-which(nas >= 3),]

naRef <- apply(calldata, 1, function(x){
  is.na(x["BFMI861-S1"])
})
map      <- map[-which(naRef),]
calldata <- calldata[-which(naRef),]

S1Ref <- t(apply(calldata, 1, function(x){
  as.numeric(x == x["BFMI861-S1"])
}))

dim(S1Ref)
dim(calldata)
dim(map)
colnames(S1Ref) <- colnames(calldata)
S1Ref[which(calldata == 1)] <- 0.5

setwd("D:/Sebastiaan/Images")
chromosomes  <- as.character(c(1:19, "X", "Y", "MT"))
cnt <- 1
for(chr in chromosomes){
  onChr <- which(map[,"Chr"] == chr)
  width = (length(onChr) * 0.03) + 100
  #width = max(100, width)
  cat(length(onChr), width,"\n")
  plotdata <- S1Ref[onChr,]
  name <- paste0("Chromosome",cnt,".png")
  if(cnt < 10) name <- paste0("Chromosome0",cnt,".png")
  png(name, width=width, height = 100)
    op <- par(las=1, mar = c(1, 5, 1, 1) + 0.1)
    image(x=1:nrow(plotdata), y=1:6, as.matrix(plotdata), xaxt='n', yaxt='n', ylab="", xlab="", col=c("yellow", "blue", "darkorange"), breaks=c(-1, 0.1, 0.6, 1), las=2)
    #axis(2,at=1:6, colnames(plotdata),las=2)
    mtext(chr, side=2, line=1, cex=2)
    box()
    abline(h = 1:6 + 0.5)
  dev.off()
  cnt <- cnt + 1
}

namez <- rev(colnames(plotdata))
namez[6] <- "BFMI861-S1 (reference)"

png("legend1.png")
plot(1:10,t='n',xaxt='n',yaxt='n',ylab="",xlab="",bty='l')
legend("topright", namez , title = "Sample ordering",cex=2)
dev.off()

png("legend1a.png")
plot(1:10,t='n',xaxt='n',yaxt='n',ylab="",xlab="",bty='l')
legend("topright", namez , title = "(top to bottom)",cex=2)
dev.off()

png("legend2.png",width=800)
plot(1:10,t='n',xaxt='n',yaxt='n',ylab="",xlab="",bty='l')
legend("topright", c("Equal to BFMI861-S1 (reference)","Different from BFMI861-S1 (reference)", "Heterozygous"), fill=c("darkorange","yellow","blue"),title = "Color legend",cex=2)
dev.off()
