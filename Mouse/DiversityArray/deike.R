# Deike.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2017
# first written Jan, 2017
#
# Analysis of Diversity array data for Deike

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")

#Deike 1
deike1  <- read.table("Analysis/Diabetes/DEIKE_BFMI861-S2vs861S1n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))

#Deike 2
deike2  <- read.table("Analysis/Diabetes/DEIKE_BFMI861-S2vs861S1_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))

#Deike 3
deike3  <- read.table("Analysis/Diabetes/DEIKE_BFMI861-S1vs861S2n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))


#Deike 3
#deike3 <- read.table(file="Analysis/Diabetes/BFMI861-S2nB6Nvs861S1n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE)


#Deike 4
#deike4  <- read.table("Analysis/Diabetes/BFMI861-S1vs861S2n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))



findRegions <- function(snps, snpInRegion = 35, maxDistance = 4000000){
  chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

  snpInRegion <- 35                                                         # A region should contain at least 35 SNPs
  maxDistance <- 4000000                                                    # If the neighbouring SNP is less that 4mb away, we continue with our region

  regions <- NULL
  for(chr in chromosomes){
    onchr  <- snps[snps[,"Chr"] == chr,]
    if(nrow(onchr) > snpInRegion){
      sorted <- sort(as.numeric(onchr[,"Location"]), index.return=TRUE)
      onchr  <- onchr[sorted$ix,]
      sloc   <- 0
      snpcnt <- 0
      for(snp in 1:(nrow(onchr)-1)){
        cloc <- onchr[snp, "Location"]
        nloc <- onchr[snp+1, "Location"]
        if(as.numeric(nloc) - as.numeric(cloc) < maxDistance){
          if(snpcnt == 0) sloc <- cloc
          snpcnt <- snpcnt + 1
        }else{
          if(snpcnt > snpInRegion){
            regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
            cat("New region:", chr, sloc, cloc, snpcnt,"\n")
          }
          snpcnt <- 0
        }
      }
      if(snpcnt > snpInRegion){
        regions <- rbind(regions, as.numeric(c(chr, sloc, cloc, snpcnt)))
        cat("New region:",chr, sloc, cloc, snpcnt,"\n")
      }
    }
  }

  colnames(regions) <- c("Chr","Start","End","SNPs")
  return(regions)
}

regionsD1 <- findRegions(deike1)
write.table(regionsD1, "Analysis/Diabetes/S2vsOther_Regions.txt", sep="\t", row.names=FALSE)

regionsD2 <- findRegions(deike2)
write.table(regionsD2, "Analysis/Diabetes/S1vsS2_Regions.txt", sep="\t", row.names=FALSE)

regionsD3 <- findRegions(deike3)
write.table(regionsD3, "Analysis/Diabetes/S1vsOther_Regions.txt", sep="\t", row.names=FALSE)

### Regions D1 = S2 vs Other
query.regions <- paste(regionsD1[,"Chr"],regionsD1[,"Start"],regionsD1[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]


write.table(res.biomart, "Analysis/Diabetes/S2vsOther_GenesInRegions.txt", sep="\t", row.names=FALSE)

### Regions D2 = S1 vs S2
query.regions <- paste(regionsD2[,"Chr"],regionsD2[,"Start"],regionsD2[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]


write.table(res.biomart, "Analysis/Diabetes/S1vsS2_GenesInRegions.txt", sep="\t", row.names=FALSE)


### Regions D3 = S1 vs Other
query.regions <- paste(regionsD3[,"Chr"],regionsD3[,"Start"],regionsD3[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]


write.table(res.biomart, "Analysis/Diabetes/S1vsOther_GenesInRegions.txt", sep="\t", row.names=FALSE)




png(file="D:/Deike3Chr.png",width=400, height= -100 + 3*400)
op <- par(mfrow=c(3,1))
op <- par(mar = c(2, 3, 3, 1) + 0.1)
op <- par(cex = 1)

#### S1

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
snpOUT       <- read.table("Analysis/Diabetes/DEIKE_BFMI861-S1vs861S2n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
regions      <- read.table("Analysis/Diabetes/S1vsOther_Regions.txt", sep="\t", header=TRUE, colClasses=c("character"))
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)

selectedmarkers <- unlist(lapply(regions[,"markersInRegions"], strsplit,", "))

mlength <- max(chrInfo[,"Length"])

#

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="BFMI861-S1 versus the other BFMI", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(snpOUT, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  #cat(x["Chr"],"->",yloc,"\n")
  points(x=xloc, y=yloc+0.25, pch='▼', col='red',cex=0.5)
})

aa <- apply(markers, 1,function(x){
  if(x["markerID"] %in% selectedmarkers){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc-0.15, pch='▲', col='blue',cex=0.5)
  }
})

aa <- apply(regions, 1, function(x){
  yloc <- match(x["Chr"], chromosomes); xlocS <- x["Start"]; xlocE <- x["End"]
  lines(c(as.numeric(xlocS)-2500000,as.numeric(xlocE)+2500000), c(yloc, yloc), type="l", col="black", lty=1,lwd=2)
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 15000000)/1000000, at=seq(0, mlength, 15000000), cex.axis=0.7)


#### S2


chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
snpOUT       <- read.table("Analysis/Diabetes/DEIKE_BFMI861-S2vs861S1n860-12n860-S2_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
regions      <- read.table("Analysis/Diabetes/S2vsOther_Regions.txt", sep="\t", header=TRUE, colClasses=c("character"))
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)

selectedmarkers <- unlist(lapply(regions[,"markersInRegions"], strsplit,", "))

mlength <- max(chrInfo[,"Length"])

#png(file="Analysis/Figures/S2vsOther_Regions.png",width=1200, height=800)

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="BFMI861-S2 versus the other BFMI", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(snpOUT, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  #cat(x["Chr"],"->",yloc,"\n")
  points(x=xloc, y=yloc+0.25, pch='▼', col='red',cex=0.5)
})

aa <- apply(markers, 1,function(x){
  if(x["markerID"] %in% selectedmarkers){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc-0.15, pch='▲', col='blue',cex=0.5)
  }
})

aa <- apply(regions, 1, function(x){
  yloc <- match(x["Chr"], chromosomes); xlocS <- x["Start"]; xlocE <- x["End"]
  lines(c(as.numeric(xlocS)-2500000,as.numeric(xlocE)+2500000), c(yloc, yloc), type="l", col="black", lty=1,lwd=2)
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 15000000)/1000000, at=seq(0, mlength, 15000000), cex.axis=0.7)




#### S1 vs S2

chromosomes  <- as.character(c(1:19, "X", "Y", "M"))

setwd("D:/Edrive/Mouse/DNA/DiversityArray/")
chrInfo      <- read.table("Annotation/mouseChrInfo.txt", header=TRUE)
snpOUT       <- read.table("Analysis/Diabetes/DEIKE_BFMI861-S2vs861S1_SNPs.txt", sep="\t", header=TRUE, colClasses=c("character"))
regions      <- read.table("Analysis/Diabetes/S2vsS1_Regions.txt", sep="\t", header=TRUE, colClasses=c("character"))
markers      <- read.table("Annotation/GeneticMarkers.txt", sep="\t", header=TRUE)

selectedmarkers <- unlist(lapply(regions[,"markersInRegions"], strsplit,", "))

mlength <- max(chrInfo[,"Length"])

#png(file="Analysis/Figures/S2vsOther_Regions.png",width=1200, height=800)

plot(c(0, mlength), c(1,nrow(chrInfo)), t='n', main="BFMI861-S1 versus BFMI861-S2", yaxt="n", ylab="Chromosome", xlab="Length (Mb)", xaxt="n")
cnt <- 1
aa <- apply(chrInfo,1,function(x){
  lines(c(0,x["Length"]), c(cnt, cnt), type="l", col="black", lty=1)
  cnt <<- cnt + 1
})

aa <- apply(snpOUT, 1,function(x){
  yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
  #cat(x["Chr"],"->",yloc,"\n")
  points(x=xloc, y=yloc+0.25, pch='▼', col='red',cex=0.5)
})

aa <- apply(markers, 1,function(x){
  if(x["markerID"] %in% selectedmarkers){
    yloc <- match(x["Chr"], chromosomes); xloc <- x["Location"]
    points(x=xloc, y=yloc-0.15, pch='▲', col='blue',cex=0.5)
  }
})

aa <- apply(regions, 1, function(x){
  yloc <- match(x["Chr"], chromosomes); xlocS <- x["Start"]; xlocE <- x["End"]
  lines(c(as.numeric(xlocS)-2500000,as.numeric(xlocE)+2500000), c(yloc, yloc), type="l", col="black", lty=1,lwd=2)
})

axis(2,chrInfo[,1], at=c(1:nrow(chrInfo)), las=1)
axis(1, seq(0, mlength, 15000000)/1000000, at=seq(0, mlength, 15000000), cex.axis=0.7)


dev.off()


regions <- findRegions(deike)
query.regions <- paste(regions[,"Chr"],regions[,"Start"],regions[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]

write.table(res.biomart, "Analysis/Diabetes/BFMI861-S2nB6Nvs861S1n860-12n860-S2_GenesInRegions.txt", sep="\t", row.names=FALSE)
write(res.biomart[,"ensembl_gene_id"], "Analysis/Diabetes/BFMI861-S2nB6Nvs861S1n860-12n860-S2_MGIsymbols.txt")




regions <- findRegions(deike1)

query.regions <- paste(regions[,"Chr"],regions[,"Start"],regions[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]

write.table(res.biomart, "Analysis/Diabetes/BFMI861-S2vsAll_GenesInRegions.txt", sep="\t", row.names=FALSE)


regions <- findRegions(deike3)

query.regions <- paste(regions[,"Chr"],regions[,"Start"],regions[,"End"], sep=":")                                  # Regions encoded for biomart

library(biomaRt)
bio.mart <- useMart(biomart="ensembl", dataset="mmusculus_gene_ensembl")                                            # Biomart for mouse genes
res.biomart <- getBM(attributes = c("ensembl_gene_id", "mgi_id","mgi_symbol", "mgi_description", "chromosome_name", "start_position", "end_position", "strand"), 
                        filters = c("chromosomal_region", "biotype"), 
                        values = list(query.regions,"protein_coding"), mart = bio.mart)

sortnames <- c("chromosome_name", "start_position")                                                                 # Sort on multiple columns
res.biomart <- res.biomart[do.call("order", res.biomart[sortnames]), ]

write.table(res.biomart, "Analysis/Diabetes/BFMI861-S1vsAll_GenesInRegions.txt", sep="\t", row.names=FALSE)


