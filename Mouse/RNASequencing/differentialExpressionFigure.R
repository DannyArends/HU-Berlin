# differentialExpressionFigure.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014
#
# Create a figure for the RNA sequencing data (pre-processed by MDC)

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

setwd("E:/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)
  
setwd("E:/Mouse/RNA/Sequencing/Reciprocal Cross B6 BFMI by MPI/")
RPKM <- read.table("Analysis/BFMI_RPKM_Qnorm_ANN_AddDom.txt", sep="\t", header=TRUE, colClasses="character")

mlength <- max(chrInfo[,"Length"])

png("MaternalOriginOurAnalysis_norm_above1RPKM.png", width = 2000, height = 1000)
  ids <<- 1
  selected <- NULL

  op <- par(cex = 2.5)
  plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="Dominant expression compared to the original strains", yaxt="n", xlab="Chromosome", ylab="Length (Mb)", xaxt="n")

  abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")

  cnt <- 1
  aa <- apply(chrInfo,1,function(x){
    lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1,lwd=3)
    cnt <<- cnt + 1
  })

  aa <- apply(RPKM, 1,function(x){
    xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
    #ttest <- as.numeric(x["tTest.1"])
    #if(!is.na(ttest) && ttest < 0.1){
      if(x["Mean.BFMI860.12xB6N.L"] > 1){
        
        col <- "white"
        if(x["A.D_BFMI860.12xB6N"] == "B6N") col <- "gray"
        if(x["A.D_BFMI860.12xB6N"] == "BFMI") col <- "orange"
        if(x["A.D_BFMI860.12xB6N"] == "ADDITIVE") col <- "white"
        if(col != "white"){
          selected <<- c(selected, ids)
          points(x=xloc + 0.15, y=yloc, pch="-", col=col,cex=1.2)
        }
      }
      if(x["Mean.B6NxBFMI860.12.L"] > 1){
        col <- "white"
        if(x["A.D_B6NxBFMI860.12"] == "B6N") col <- "gray"
        if(x["A.D_B6NxBFMI860.12"] == "BFMI") col <- "orange"
        if(x["A.D_B6NxBFMI860.12"] == "ADDITIVE") col <- "white"
        if(col != "white"){
          selected <<- c(selected, ids)
          points(x=xloc - 0.15, y=yloc, pch="-", col=col,cex=1.2)
        }
      }
      
      ids <<- ids + 1
    #}
  })

  axis(1,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
  axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7, las=2)
  legend("topright", c("BFMI like", "B6N like"), fill=c("orange","gray"))
dev.off()

cat("Maternal BFMI: B6N:", sum(RPKM[,"A.D_BFMI860.12xB6N"] == "B6N"), "BFMI:", sum(RPKM[,"A.D_BFMI860.12xB6N"] == "BFMI"),"\n")
cat("Maternal B6: B6N:", sum(RPKM[,"A.D_B6NxBFMI860.12"] == "B6N"), "BFMI:", sum(RPKM[,"A.D_B6NxBFMI860.12"] == "BFMI"),"\n")

png("OriginOfExpressionOurAnalysis_norm_above1RPKM.png", width = 2000, height = 1000)
  ids <<- 1
  maternal <- NULL
  paternal <- NULL
  op <- par(cex = 2.5)
  plot(y=c(0, mlength), x=c(1,nrow(chrInfo)), t='n', main="Origin of Expression", yaxt="n", ylab="Length (Mb)", xlab="Chromosome", xaxt="n")

  abline(h=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")
  cnt <- 1
  aa <- apply(chrInfo,1,function(x){
    lines(c(cnt, cnt), c(0,x["Length"]), type="l", col="gray", lty=1,lwd=3)
    cnt <<- cnt + 1
  })

  aa <- apply(RPKM, 1,function(x){
    if(x["Mean.BFMI860.12xB6N.L"] > 1 && x["Mean.B6NxBFMI860.12.L"] > 1){
      xloc <- match(as.character(x["chromosome_name"]), chromosomes); yloc <- as.numeric(x["start_position"])
      if(x["A.D_BFMI860.12xB6N"] == "B6N" && x["A.D_B6NxBFMI860.12"] == "BFMI"){
        points(x=xloc, y=yloc, pch="-", col="blue",cex=1.5)
        paternal <<- c(paternal, ids)
      }
    
      if(x["A.D_BFMI860.12xB6N"] == "BFMI" && x["A.D_B6NxBFMI860.12"] == "B6N"){
        points(x=xloc, y=yloc, pch="-", col="red",cex=1.5)
        maternal <<- c(maternal, ids)
      }
    }
    ids <<- ids + 1
  })

  axis(1,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
  axis(2, seq(0, mlength, 10000000)/1000000, at=seq(0, mlength, 10000000), cex.axis=0.7, las=2)
  legend("topright", c("Expressed as paternal", "Expressed as maternal"), fill=c("blue","red"))
  # TODO: make a supplementary table of these genes.
dev.off()

RPKM <- cbind(RPKM, "DominanceEffects" = rep(NA, nrow(RPKM)))
RPKM[maternal,"DominanceEffects"] <- "Maternal"
RPKM[paternal,"DominanceEffects"] <- "Paternal"

smallRPKM <- RPKM[unique(selected),]
ids <- 1
aa <- apply(smallRPKM, 1, function(x){
  if(x["A.D_BFMI860.12xB6N"] == "BFMI" && x["A.D_B6NxBFMI860.12"] == "BFMI"){
    smallRPKM[ids, "DominanceEffects"] <<- "BFMI"
  }else if(x["A.D_BFMI860.12xB6N"] == "B6N"  && x["A.D_B6NxBFMI860.12"] == "B6N"){
    smallRPKM[ids, "DominanceEffects"] <<- "B6N"
  }
  ids <<- ids +1
})
smallRPKM[is.na(smallRPKM[,"DominanceEffects"]),"DominanceEffects"] <- "Complex"

write.table(smallRPKM, "RPKM_Filtered_Dominance.txt", sep="\t", row.names=FALSE)