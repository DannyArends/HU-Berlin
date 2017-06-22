

source("D:/Ddrive/Github/HU-Berlin/Mouse/Muga/ATB_Paper/dateToSeason.R")

setwd("D:/Edrive/Mouse/DNA/MegaMuga/")                                                                               # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE, colClasses="character")      # Normal A, H, B genotypes
phenotypes <- read.csv("Phenotypes/MatchedPhenotypes.txt", sep="\t", header=TRUE)

F2 <- c(rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)], "BFMI860-12 (V2)")                                                # The F2 individuals

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"]))                                     # Add the season column to the matrix
birthmonth <- unlist(lapply(strsplit(as.character(phenotypes[,"W.dat"]),".", fixed=TRUE),"[",2))
phenotypes <- cbind(phenotypes, Birthmonth = birthmonth)                                                      # Add the birth month column to the matrix

phenotypes <- cbind(phenotypes, mri42d_fatDlean = phenotypes[,"mri42d_fat"] / phenotypes[,"mri42d_lean"])     # Fat / Lean day 42
phenotypes <- cbind(phenotypes, mri56d_fatDlean = phenotypes[,"mri56d_fat"] / phenotypes[,"mri56d_lean"])     # Fat / Lean day 56
phenotypes <- cbind(phenotypes, mri70d_fatDlean = phenotypes[,"mri70d_fat"] / phenotypes[,"mri70d_lean"])     # Fat / Lean day 70



from <- which(rownames(numericgenotypes) == "UNC5009661")
to <- which(rownames(numericgenotypes) == "JAX00519845") + 5


#fatMice <- rownames(phenotypes[which(phenotypes[,"mri70d_fat"] > 12),])
#leanMice <- rownames(phenotypes[which(phenotypes[,"mri70d_fat"] < 7),])

#noBFMIinShoulder <- colnames(genotypes[,which(genotypes["UNC5035770",] != "B")])

#genotypes <- genotypes[,noBFMIinShoulder]
#phenotypes <- phenotypes[noBFMIinShoulder,]

#qtls <- mriGWAS(genotypes,   phenotypes, "mri42d_fatDlean");                                                                # Map QTLs, normal GWAS (A, H, B)

#from <- which(rownames(qtls) == "UNC5009661")
#to <- which(rownames(qtls) == "JAX00519845")

#plot(qtls[from:to,"marker"], t = 'l')

BFMIish <- apply(genotypes[from:to,], 2,function(x){
  x == genotypes[from:to,"BFMI860-12 (V2)"]
})

BFMIish <- apply(BFMIish,2,as.numeric)
BFMIish[genotypes[from:to,] == "H"] <- 0.5

aaa <- c("black", "gray", "orange")
names(aaa) <- c(1,3,5)

locations <- as.numeric(map[rownames(numericgenotypes[from:to,]),"Mb_NCBI38"])

locs <- NULL
mstart <- locations[1]
for(x in 2:length(locations)){
  if(mstart != locations[x-1]){ mstart = mstop; cat(".") }
  mstop <- mean(c(locations[x-1], locations[x]))
  if(x == length(locations)) mstop <- locations[x]
  cat(mstart, mstop,"\n")
  locs <- rbind(locs, c(mstart, mstop))
}

op <- par(mar=c(5, 4, 4, 4) + 0.1)
fat <- rownames(phenotypes[which(phenotypes[,"mri70d_fatDlean"] > 0.5),])
lean <- rownames(phenotypes[which(phenotypes[,"mri70d_fatDlean"] < 0.1),])
N <- length(fat) + length(lean)
plot(c(34111648, max(as.numeric(map[rownames(numericgenotypes[from:to,]),"Mb_NCBI38"]))), c(1,N+1), t = 'n', xlab="Distance in base pairs",yaxt='n',ylab="Weight")

cnt <- 1
for(ind in lean){
  points(x=locations, y = rep(cnt,length(from:to)), col=aaa[1 + 2*numericgenotypes[from:to,ind]],pch=15)
  for(x in 1:nrow(locs)){ lines(c(locs[x,1],locs[x,2]),c(cnt,cnt), col=aaa[1 + 2*numericgenotypes[from:to,ind][x]]); }
  cnt <- cnt + 1
}

#points(x=locations, y = rep(2,length(from:to)), col=aaa[1 + 2*numericgenotypes[from:to,"6661155"]],pch=15)
#for(x in 1:nrow(locs)){ lines(c(locs[x,1],locs[x,2]),c(2,2), col=aaa[1 + 2*numericgenotypes[from:to,"6661155"][x]]); }
#points(x=locations, y = rep(3,length(from:to)), col=aaa[1 + 2*numericgenotypes[from:to,"6661907"]],pch=15)
#for(x in 1:nrow(locs)){ lines(c(locs[x,1],locs[x,2]),c(3,3), col=aaa[1 + 2*numericgenotypes[from:to,"6661907"][x]]); }
cnt <- cnt + 1

for(ind in fat){
  points(x=locations, y = rep(cnt,length(from:to)), col=aaa[1 + 2*numericgenotypes[from:to,ind]],pch=15)
  for(x in 1:nrow(locs)){ lines(c(locs[x,1],locs[x,2]),c(cnt,cnt), col=aaa[1 + 2*numericgenotypes[from:to,ind][x]]); }
  cnt <- cnt + 1
}
  


#points(x=locations, y = rep(5,length(from:to)), col=aaa[1 + 2*numericgenotypes[from:to,"6661817"]],pch=15)
#for(x in 1:nrow(locs)){ lines(c(locs[x,1],locs[x,2]),c(5,5), col=aaa[1 + 2*numericgenotypes[from:to,"6661817"][x]]); }
#points(x=locations, y = rep(6,length(from:to)), col=aaa[1 + 2*numericgenotypes[from:to,"6661154"]],pch=15)
#for(x in 1:nrow(locs)){ lines(c(locs[x,1],locs[x,2]),c(6,6), col=aaa[1 + 2*numericgenotypes[from:to,"6661154"][x]]); }
#points(x=locations, y = rep(7,length(from:to)), col=aaa[1 + 2*numericgenotypes[from:to,"6661872"]],pch=15)
#for(x in 1:nrow(locs)){ lines(c(locs[x,1],locs[x,2]),c(7,7), col=aaa[1 + 2*numericgenotypes[from:to,"6661872"][x]]); }

abline(v=36808838, lty=2)
abline(v=36537006, lty=2)

axis(2, 1:(N+1), c(rep("lean", length(lean)),"",rep("fat", length(fat))),las=2,cex.axis=0.6)
axis(4, 1:(N+1), c(lean,"",fat),las=2,cex.axis=0.6)



