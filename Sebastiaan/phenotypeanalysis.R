# Spider plot, create the axis out of different traits
# Use a circle to align, every axis is scaled from: min = 0, max = sqrt(2), and from [0,0] plotted at angle rot
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Sep, 2015
# first written Sep, 2015

setwd("D:/Collegues/Sebastiaan")
phenotypes <- read.table("All_Measurements_20weeks.txt", sep = "\t", header = TRUE, na.strings = c(NA, ".", "/"))
phenames <- colnames(phenotypes)[c(6:ncol(phenotypes))]                                           # Take the ones we can plot, Malonyl_CoA can't be sued

clusters <- hclust(dist(cor(phenotypes[,phenames], use = "pair", method = "spearman")))           # Cluster phenotypes based on correlation similarity
#phenames <- phenames[clusters$order]                                                              # Order them in a 'logical' way

phenames <- c("FAT140", "FATpro140","LEAN140", "IntramuscFat_Q", 
              "IntramuscFat_LD", "Liver_TRIGS", "Insulin","ITT20_0", "ITT20_15","ITT20_30", "ITT20_60", "AUC20", "C", "Akt1", "Cbl", 
              "Foxa2" ,"Irs1", "Irs2", "Insr", "Igf1r","Lep", "Prkaa2", "Slc2a4")

mline <- as.character(phenotypes[,"Line"])
mline[which(mline == "852_S1")] <- "852"
mline[which(mline == "856_S2")] <- "856"

phenotypes[,"Line"] <- mline


rot <- seq(1, 360, 360/length(phenames))                                                          # Rotational axis
names(rot) <- phenames                                                                            # Phenotype names displayed for each axis
cp <- seq(0, 2 * pi, length = 360)                                                                # Circle points in degrees
coords <- t(rbind( sin(cp), cos(cp)))                                                             # Circle coordinates

sublines <- c("852", "856", "860_12", "860_S2", "861_S1", "861_S2")
refs <- c("B6_J", "DBA_J")
colz <- c(rgb(1,0,0,0.2), rgb(0,1,0,0.2), rgb(0,0,1,0.2))                                         # B6 = Red,  DBA = Green, Subline = Blue
colzNT <- c(rgb(1,0,0), rgb(0,1,0), rgb(0,0,1))                                                   # B6 = Red,  DBA = Green, Subline = Blue


tiff("Heise_Fig5.tif", width = 480 * 5, height = 480 * 6, res=300, compression = "lzw")
op <- par(mfrow=c(3, 2), mai = c(0.05, 0.1, 0.2, 0.1))                                            # 6 sublines ( so we use a 3x2 plot )

for(y in 1:length(sublines)){
  plot(c(-1.2, 1.2),c(-1.2, 1.2), t = 'n', ylab="", xlab="",xaxt='n',yaxt='n', main=paste0("BFMI",gsub("_", "-", sublines[y])))  # Create the plot window
  points(coords, pch = 18, cex = 0.4)                                                             # Add the circle
  for(x in 1:length(rot)){
    points(c(0, coords[rot[x],1]), c(0,coords[rot[x],2]), t='l', lwd=0.5, lty=3)                  # The different phenotype axis
    text(1.1 * coords[rot[x],1], 1.1 * coords[rot[x],2], names(rot)[x], cex=0.7)                  # The name on the outside of the circle
  }
  cnt <- 1
  for(subline in c(refs, sublines[y])){                                                           # 2 reference lines + the selected line
    scoords <- NULL
    for(phe in names(rot)){
      mL    <- mean(phenotypes[phenotypes[,"Line"] == subline, phe], na.rm=TRUE)                  # Mean phenotype for this 'subline'
      minL  <- min(phenotypes[, phe],na.rm=TRUE) ; maxL  <- max(phenotypes[, phe],na.rm=TRUE)     # Min and max across all
      
      if(is.nan(mL)) stop(paste("Not enough observations for trait:", phe))                       # Warn if we cannot get a mean (we have to remove the trait then)
      phescale <- ((mL - minL) / (maxL - minL))                                                   # Scale
      coord <- c(phescale * coords[rot[phe],1], phescale * coords[rot[phe],2])                    # X and Y coordinate on the scaled circle
      points(t(coord), pch=19,cex=0.8, col = colzNT[cnt])                                         # plot the dot
      scoords <- rbind(scoords, coord)                                                            # remember the coordinate
    }
    rownames(scoords) <- names(rot)                                                               # Set rownames on the polygon edges
    polygon(scoords, col = colz[cnt], pch=18,fillOddEven = TRUE)                                  # Connect all the dots using a polygon
    cat("Done with subline", cnt,"\n")
    cnt <- cnt + 1
  }
}
dev.off()
