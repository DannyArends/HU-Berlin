# Analysis of JAX data from the mouse diversity array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014
#

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")
library(ape)

calldata <- read.table(file="Analysis/measurementsALL_annotated.txt", sep="\t", header=TRUE)

clusters   <- hclust(dist(t(calldata[,9:ncol(calldata)])))
dendrogram <- as.dendrogram(clusters)                                                             # Dendrogram of all ALL lines
phylogram  <- as.phylo(clusters)                                                                  # Phylogram of all ALL lines

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                           # Load the annotation
annInData  <- match(colnames(calldata), rownames(annotation))                                     # Align JAX B6 with the JAX BMMI data
annotation <- annotation[annInData,]

strain2color <- cbind(levels(annotation[,"Strain"]), rainbow(length(levels(annotation[,"Strain"]))))
colnames(strain2color) <- c("Strain", "Color")
rownames(strain2color) <- strain2color[,1]

# Function to color the dendrogram nodes based on the strain
labelCol <- function(x){
  if(is.leaf(x)){
    label <- attr(x, "label")                                            # Fetch label
    strain <- annotation[label,"Strain"]                                 # Convert the label to the strain
    attr(x, "nodePar") <- list(lab.col = strain2color[strain, "Color"])  # Set the label color based on the strain
  }
  return(x)
}

dendroColor <- dendrapply(dendrogram, labelCol)
png(file="Analysis/Figures/DendrogramALL-lines.png",width=1200, height=800)
  par(cex=0.8, font=1, mai=c(1.5,1,1,1))                                                          # Set some parameters for the plot
  plot(dendroColor, type="triangle", main="Dendrogram ALL lines", xlab="", center=TRUE)
  legend("topleft",rownames(strain2color), col=strain2color[,"Color"], lwd=2, cex=0.8)
dev.off()

png(file="Analysis/Figures/PhyloFanALL-lines.png",width=800, height=800)
  plot(phylogram, type = "fan", tip.col=strain2color[annotation[phylogram$tip.label,"Strain"],"Color"])
  legend("topleft",rownames(strain2color), col=strain2color[,"Color"], lwd=2, cex=0.8)
dev.off()



s1 <- "X8614310"
s2 <- "X8614277"

sum(apply(calldata[, c(s1,s2)],1,function(x){
  x[1]!=x[2]
}),na.rm=TRUE)


