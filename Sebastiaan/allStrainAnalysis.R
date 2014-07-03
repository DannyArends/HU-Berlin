# Analysis of JAX data from the mouse diversity array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Juli, 2014
# first written Juli, 2014
#

setwd("E:/Mouse/DiversityArray/")

calldata <- read.table(file="Analysis/measurementsALL_annotated.txt", sep="\t", header=TRUE)

dendrogram <- as.dendrogram(hclust(dist(t(calldata[,9:ncol(calldata)]))))                       # Dendrogram of all ALL lines

annotation <- read.table("Annotation/MouseAnnotation.txt", header=TRUE)                           # Load the annotation
annInData  <- match(colnames(calldata), rownames(annotation))                                    # Align JAX B6 with the JAX BMMI data
annotation <- annotation[annInData,]

strain2color <- cbind(levels(annotation[,"Strain"]), rainbow(length(levels(annotation[,"Strain"]))))
colnames(strain2color) <- c("Strain", "Color")
rownames(strain2color) <- strain2color[,1]

# Function to set label color to strains
labelCol <- function(x) {
  if(is.leaf(x)){
    label <- attr(x, "label")                                            # Fetch label
    strain <- annotation[label,"Strain"]                                 # Fetch Strain
    attr(x, "nodePar") <- list(lab.col = strain2color[strain, "Color"])  # Set label color based on strain
  }
  return(x)
}

dendroColor <- dendrapply(dendrogram, labelCol)
png(file="Analysis/Figures/DendrogramALL-lines.png",width=1200, height=800)
  par(cex=0.8,font=1, mai=c(1.5,1,1,1))
  plot(dendroColor, main="Dendrogram ALL lines", xlab="", center=TRUE)
  legend("topleft",rownames(strain2color), col=strain2color[,"Color"], lwd=2, cex=0.8)
dev.off()
