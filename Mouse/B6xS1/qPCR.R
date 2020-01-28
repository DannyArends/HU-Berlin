#
# RT pcr analysis Sandra
#

setwd("D:/Edrive/Mouse/AIL_Manuel")
mdata <- read.csv("qPCR_Sandra.txt", sep = "\t", na.strings = c("Undetermined", ""), colClasses = "character")
# Remove things with missing CT values
mdata <- mdata[-which(is.na(mdata[, "CT"])),]

# Fix some issues
mdata[,"Target.Name"] <- tolower(mdata[,"Target.Name"])
mdata[which(mdata[, "Sample.Name"] %in% c("primt0,1", "primt1", "primt10")),"Sample.Name"] <- "primertest"

# Our marker / Genotype data
marker <- "S1H083826428"
marker.data <- read.table("genomatrix.clean.txt", sep = "\t", colClasses = "character")
genotypes <- marker.data[marker,]
names(genotypes) <- gsub("X666.", "", names(genotypes), fixed=TRUE)
genotypes <- genotypes[which(names(genotypes) %in% mdata[, "Sample.Name"])]

### Group into genotypes by the marker near the CES genes
mdata <- mdata[which(mdata[, "Sample.Name"] %in% names(genotypes)),]

mdata <- cbind(mdata, GT = as.character(genotypes[mdata[, "Sample.Name"]]))

samples <- unique(mdata[, "Sample.Name"])
genes <- unique(mdata[, "Target.Name"])
housekeeper <- c("actb", "rsp")
genes <- genes[which(!(genes %in% housekeeper))]

### Copy the original data, so we can compare back
mdataOld <- mdata

redo <- c() # samples that need to be redone
# Fix/Remove the large CT.SD values
for (g in c(genes,housekeeper)) {
  for (s in samples) {
    idx <- which(mdata[, "Sample.Name"] == s &  mdata[, "Target.Name"] == g)
    gData <- mdata[idx,]
    vals <- na.omit(as.numeric(gData[, "CT"]))
    CTmean <- round(mean(vals, na.rm = TRUE), 1)
    CTsd <- round(sd(vals, na.rm = TRUE), 1)
    if (is.na(CTsd) || CTsd > 0.2) {
      cat(s, " ", g, " CTSD is too large: ", CTsd, ":", as.numeric(gData[, "CT"]), ", ")
      # If we have 3 values figure out if removing the one furthest from the mean makes the data suitable
      if (length(vals) == 3) {
        biggestDif <- which.max(abs(vals - CTmean))
        CTsdC <- round(sd(vals[-biggestDif], na.rm = TRUE), 1)
        if (CTsdC <= 0.2) {
          cat("Fix one, removed:",vals[biggestDif],"\n")
          mdata[idx[biggestDif],"CT"] <- NA
        } else {
          cat("Removing ALL\n")
          redo <- rbind(redo, c(s,g))
          mdata[idx,"CT"] <- NA
        }
      } else {
        cat("No 3 measurements, Removing ALL\n")
        redo <- rbind(redo, c(s,g))
        mdata[idx,"CT"] <- NA
      }
    }
  }
}
write.table(redo, "failed_qPCR_toREDO.txt", sep = "\t", quote = FALSE, row.names=FALSE)

# Compute CT values relative to the housekeepers
pvals <- c()
for(g in genes){
  mymatrix <- c()
  for(s in samples){
    gData <- mdata[which(mdata[, "Sample.Name"] == s &  mdata[, "Target.Name"] == g),]
    gMean <- round(mean(as.numeric(gData[, "CT"]), na.rm=TRUE),1)
    myrow <- c(g, genotypes[s], gMean)
    for(hk in housekeeper){
      hKeeper <- mdata[which(mdata[, "Sample.Name"] == s &  mdata[, "Target.Name"] == hk),]
      hKeeperMean <- round(mean(as.numeric(hKeeper[, "CT"]), na.rm=TRUE),1)
      myrow <- c(myrow, hKeeperMean)
    }
    mymatrix <- rbind(mymatrix, myrow)
  }
  rownames(mymatrix) <- samples
  colnames(mymatrix) <- c("Gene", "GT", "CT", housekeeper)
  dCTlist <- vector("list", length(unique(mymatrix[, "GT"])))
  names(dCTlist) <- unique(mymatrix[, "GT"])
  for(gt in unique(mymatrix[, "GT"])){
    gtmatrix <- mymatrix[which(mymatrix[, "GT"] == gt),]
    CTdiff_actb <- as.numeric(gtmatrix[, "CT"]) - as.numeric(gtmatrix[, "actb"])
    CTdiff_actb <- 2^(-CTdiff_actb)
    CTdiff_rsp <- as.numeric(gtmatrix[, "CT"]) - as.numeric(gtmatrix[, "rsp"])
    CTdiff_rsp <- 2^(-CTdiff_rsp)
    mm <- cbind(CTdiff_actb, CTdiff_rsp)
    rownames(mm) <- rownames(gtmatrix)
    dCTlist[[gt]] <- mm
  }
  for(gtA in c("AA", "AG")){
    for(gtB in c("AG", "GG")){
      if(gtA != gtB){
        P_actb <- t.test(dCTlist[[gtA]][, "CTdiff_actb"], dCTlist[[gtB]][, "CTdiff_actb"])$p.value
        P_rsp <- t.test(dCTlist[[gtA]][, "CTdiff_rsp"], dCTlist[[gtB]][, "CTdiff_rsp"])$p.value
        pvals <- rbind(pvals, c(g, gtA, gtB, "actb", round(mean(dCTlist[[gtA]][, "CTdiff_actb"],na.rm=TRUE),3), round(sd(dCTlist[[gtA]][, "CTdiff_actb"],na.rm=TRUE),3), round(mean(dCTlist[[gtB]][, "CTdiff_actb"],na.rm=TRUE),3), round(sd(dCTlist[[gtB]][, "CTdiff_actb"],na.rm=TRUE),3), P_actb))
        pvals <- rbind(pvals, c(g, gtA, gtB, "rsp", round(mean(dCTlist[[gtA]][, "CTdiff_rsp"],na.rm=TRUE),3), round(sd(dCTlist[[gtA]][, "CTdiff_rsp"],na.rm=TRUE),3), round(mean(dCTlist[[gtB]][, "CTdiff_rsp"],na.rm=TRUE),3), round(sd(dCTlist[[gtB]][, "CTdiff_rsp"],na.rm=TRUE),3), P_rsp))
      }
    }
  }
}
colnames(pvals) <- c("Gene", "GT1", "GT2", "Housekeeper", "Mean[GT1]", "SD[GT1]", "Mean[GT2]", "SD[GT2]", "p.value")
data.frame(pvals)

colnames(redo) <- c("Sample", "Gene")
redo