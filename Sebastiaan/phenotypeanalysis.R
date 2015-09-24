setwd("D:/Sebastiaan")

phenotypes <- read.table("Tuerkei_all_measurements_20weeks.txt", sep = "\t",header = TRUE,na.strings=c(NA,"."))


rot <- c(1,135,225)
names(rot) <- c("FATpro140", "C", "IRS2")

t <- seq(0, 2 * pi,length = 360)
coords <- t(rbind( sin(t), cos(t)))

makeTransparent<-function(someColor, alpha=30)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

sublines <- c("852_S1", "856_S2", "860_12", "860_S2", "861_S1", "861_S2")
refs <- c("B6_J", "DBA_J")

op <- par(mfrow=c(3, 2))

colz <- c(rgb(1,0,0,0.1),rgb(0,1,0,0.1),rgb(0,0,1,0.1))

for(y in 1:length(sublines)){
  plot(c(-1.2, 1.2),c(-1.2, 1.2), t = 'n', ylab="", xlab="",xaxt='n',yaxt='n', main=sublines[y])
  points(coords, pch = 18,cex = 0.4)
  for(x in 1:length(rot)){
    points(c(0, coords[rot[x],1]), c(0,coords[rot[x],2]), t='l',lwd=0.5)
    text(1.1 * coords[rot[x],1], 1.1 * coords[rot[x],2], names(rot)[x])
  }
  cnt <- 1
  for(subline in c(refs, sublines[y])){
    scoords <- NULL
    for(phe in names(rot)){
      mL    <- mean(phenotypes[phenotypes[,"Line"] == subline, phe],na.rm=TRUE)
      maxL  <- max(phenotypes[, phe],na.rm=TRUE)
      minL  <- min(phenotypes[, phe],na.rm=TRUE)
      scoords <- rbind(scoords, c(((mL- minL)/(maxL - minL)) * coords[rot[phe],1], ((mL- minL)/(maxL-minL)) * coords[rot[phe],2]))
    }
    scoords <- rbind(scoords,scoords[1,])
    polygon(scoords, col = colz[cnt],pch=18)
    cat(cnt,"\n")
    cnt <- cnt + 1
  }
}


max(phenotypes[,"FATpro140"],na.rm=TRUE)