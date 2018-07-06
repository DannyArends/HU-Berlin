# Preprocessing of the MegaMuga data, mapping QTLs on different genetic maps
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

getSeason <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES),".", fixed=TRUE),"[",2)))
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5] <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8] <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11] <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2] <- "Winter"
  return(ret)
}
toAdd <- function(x){
  x[x=="A"] <- 0; x[x=="H"] <- 1; x[x=="B"] <- 2
  return(as.numeric(x))
}
toDomDev <- function(x){
  x[x=="A"] <- 0; x[x=="H"] <- 1; x[x=="B"] <- 0
  return(as.numeric(x))
}

library(parallel)

setwd("D:/Edrive/Mouse/DNA/MegaMuga/") # Read in the data from the Mega Muga
map <- read.table("Analysis/map.txt", sep="\t", colClasses=c("character"))
genotypes   <- read.table("Analysis/genotypes.txt", sep="\t", check.names=FALSE, colClasses="character") # Normal A, H, B genotypes

parentgeno <- genotypes[,c("BFMI860-12 (V2)", "B6N")]

setwd("D:/Edrive/Mouse/ClassicalPhenotypes/Reciprocal Cross B6 BFMI") # Read in the phenotypes
phenotypedata <- read.csv("20140801_AIL1_666.txt", sep="\t", header=TRUE)

phenos <- c("Vater", "W.dat", "W.Label", "Eltern_ID", "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71", "GF1", "GF2", "total.GF", "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",
            "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4", "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4", "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4", "WG", "WG2", "Farbe", "sex", "Gen.")
phenotypes <- phenotypedata[which(rownames(phenotypedata) %in% colnames(genotypes)),phenos] # Use only the phenotypes for which we have genotypes
F2 <- rownames(phenotypes)[which(phenotypes[, "Gen."] == 28)] # The F2 individuals
F2 <- F2[-which(F2=="6661459")]     

phenotypes <- cbind(phenotypes, Season = getSeason(phenotypes[,"W.dat"])) # Add the season column to the matrix
m <- "JAX00191930"

genotypes <- genotypes[,F2]
phenotypes <- phenotypes[F2,]

for(m in 1:nrow(genotypes)){
  mtable <- table(as.character(genotypes[m,]))
  small <- which(mtable < 15)
  if(length(small) > 0) genotypes[m, which(genotypes[m,] %in% names(small))] <- NA
  if(m %% 1000 == 0)cat(m, "/", nrow(genotypes), "\n")
}

onegenotype <- which(lapply(apply(genotypes, 1, table), length) <= 1) # Markers with only one genotype cannot be used in QTL mapping
genotypes   <- unique(genotypes[-onegenotype,]) # Only take the F2 individuals
cat("Left with", nrow(genotypes), "markers\n") # == Left with 10949 markers

parentgeno <- parentgeno[rownames(genotypes),]
map <- map[rownames(genotypes),]
map <- cbind(map, sumPos = NA)
chrs <- 1:20
names(chrs) <- c(1:19, "X")

chrs.starts <- c(0)
chrs.lengths <- c()
chrs.summed <- 0
chr.gap <- 25000000

for(chr in names(chrs)){
  onChr <- which(map[,"Chr"] == chr)
  chr.length <- max(as.numeric(map[onChr, "Mb_NCBI38"]))
  map[onChr,"sumPos"] <- as.numeric(map[onChr,"Mb_NCBI38"]) + chrs.summed
  chrs.summed <- chrs.summed + chr.length + chr.gap
  chrs.lengths <- c(chrs.lengths, chr.length + chr.gap)
  chrs.starts <- c(chrs.starts, chrs.summed)
}
chrs.lengths <- c(chrs.lengths, NA)

write.table(genotypes, file="cleaned_genotypes_F2.txt", sep='\t', quote=FALSE)
write.table(map, file="cleaned_map.txt", sep='\t', quote=FALSE)
write.table(phenotypes, file="cleaned_phenotypes_F2.txt", sep='\t', quote=FALSE)

# Group sizes
gttables <- apply(genotypes, 1, table)

nA <- unlist(lapply(gttables,"[", "A"))
nH <- unlist(lapply(gttables,"[", "H"))
nB <- unlist(lapply(gttables,"[", "B"))
names(nA) <- names(nH) <- names(nB) <- rownames(genotypes)

# Extract covariates
timepoints <- c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70")
subfamily <- as.factor(phenotypes[, "Vater"]) # Fixed effect: Subfamily structure (factor)
littersize <- as.numeric(phenotypes[, "WG2"]) # Fixed effect: Size of the litter (linear effect)
litternumber <- as.numeric(phenotypes[, "W.Label"] != "A") + 1
litternumber <- as.factor(litternumber) # Fixed effect: Number of litter (factor)
season <- as.factor(phenotypes[, "Season"]) # Fixed effect: Season when born (factor)
topmarker <- as.factor(unlist(genotypes["UNC5048297", ])) # topmarker 

### Map everything using a per timepoint linear model
results.lm <- vector("list", length(timepoints))
names(results.lm) <- timepoints
for(tp in timepoints){
  cl <- makeCluster(getOption("cl.cores", 6))
  clusterExport(cl, "subfamily") # export subfamily to the nodes
  clusterExport(cl, "littersize") # export littersize to the nodes
  clusterExport(cl, "litternumber") # export litternumber to the nodes
  clusterExport(cl, "season") # export season to the nodes
  clusterExport(cl, "toAdd") # export toAdd to the nodes
  clusterExport(cl, "toDomDev") # export toDomDev to the nodes
  res <- parApply(cl, genotypes, 1, function(x, phenotype){
    MA <- toAdd(as.character(x))
    MD <- toDomDev(as.character(x))
    MF <- as.factor(as.character(x))
    model.lm <- lm(phenotype ~ subfamily + littersize + litternumber + season + MA + MD)
    model.flm <- lm(phenotype ~ subfamily + littersize + litternumber + season + MF)
    pvals <- rep(NA, 7)
    names(pvals) <- c("subfamily", "littersize", "litternumber", "season", "MA", "MD", "MF")
    pvals[1:6] <- anova(model.lm)[names(pvals)[1:6],"Pr(>F)"]
    pvals[7] <- anova(model.flm)[names(pvals)[7],"Pr(>F)"]
    return(pvals)
  }, phenotype = as.numeric(phenotypes[, tp]))
  stopCluster(cl)
  rownames(res) <- c("subfamily", "littersize", "litternumber", "season", "MA", "MD","MF")
  results.lm[[tp]] <- res
  cat("Done ", tp, "\n")
}

plot(c(0, max(map[onChr, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n')
for (tp in timepoints) {
  for (chr in names(chrs)) {
    onChr <- which(map[,"Chr"] == chr)
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][5, onChr]), t ='l', col='gray',lwd=3)
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][6, onChr]), t ='l', col='blue',lwd=2)
    points(map[onChr, "sumPos"], -log10(results.lm[[tp]][7, onChr]), t ='l', col='green', lwd=1)
  }
}
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")

### Map everything using a per timepoint linear mixed model, using subfamily as a random effect
results.lmm <- vector("list", length(timepoints))
names(results.lmm) <- timepoints
for(tp in timepoints){
  cl <- makeCluster(getOption("cl.cores", 6))
  clusterExport(cl, "subfamily") # export subfamily to the nodes
  clusterExport(cl, "littersize") # export littersize to the nodes
  clusterExport(cl, "litternumber") # export litternumber to the nodes
  clusterExport(cl, "season") # export season to the nodes
  clusterExport(cl, "toAdd") # export toAdd to the nodes
  clusterExport(cl, "toDomDev") # export toDomDev to the nodes
  clusterEvalQ(cl, library(lme4)) # load lme4 per node
  res <- parApply(cl, genotypes, 1, function(x, phenotype){
    noG <- which(is.na(x))
    MA <- toAdd(as.character(x))
    MD <- toDomDev(as.character(x))
    if(length(noG) > 0){
      phenotype <- phenotype[-noG]
      MA <- MA[-noG]
      MD <- MD[-noG]
      subfamily <- subfamily[-noG]
      littersize <- littersize[-noG]
      litternumber <- litternumber[-noG]
      season <- season[-noG]
    }
    model.full <- lmer(phenotype ~ (1|subfamily) + littersize + litternumber + season + MA + MD)
    model.add <- lmer(phenotype ~ (1|subfamily) + littersize + litternumber + season + MA)
    model.dom <- lmer(phenotype ~ (1|subfamily) + littersize + litternumber + season + MD)
    model.null <- lmer(phenotype ~ (1|subfamily) + littersize + litternumber + season)
    pvals <- rep(NA, 3)
    names(pvals) <- c("Add", "Dom", "Full")
    pvals["Add"] <- anova(model.full, model.dom)["model.full", "Pr(>Chisq)"]
    pvals["Dom"] <- anova(model.full, model.add)["model.full", "Pr(>Chisq)"]
    pvals["Full"] <- anova(model.full, model.null)["model.full", "Pr(>Chisq)"]
    return(pvals)
  }, phenotype = as.numeric(phenotypes[, tp]))
  stopCluster(cl)
  rownames(res) <- c("Add", "Dom", "Full")
  results.lmm[[tp]] <- res
  cat("Done ", tp, "\n")
}

tpcolors <- colorRampPalette(c("lightblue", "darkblue"))(length(timepoints))
names(tpcolors) <- timepoints

plot(c(0, max(map[, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome", xaxt='n')
for (tp in timepoints) {
  for (chr in names(chrs)) {
    onChr <- which(map[,"Chr"] == chr)
    points(map[onChr, "sumPos"], -log10(results.lmm[[tp]]["Full",onChr]), t ='l', col = tpcolors[tp])
  }
}
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")

# Group everything into 1 big model
phenotype <- NULL
individual <- NULL
subfamily <- NULL
littersize <- NULL
litternumber <- NULL
season <- NULL
timepoint <- NULL
for(tp in timepoints){
  phenotype <- c(phenotype, phenotypes[, tp])
  individual <- c(individual, rownames(phenotypes))
  subfamily <- c(subfamily, phenotypes[, "Vater"])
  littersize <- c(littersize, phenotypes[, "WG2"])
  litternumber <- c(litternumber, as.character(phenotypes[, "W.Label"]))
  season <- c(season, as.character(phenotypes[, "Season"]))
  timepoint <- c(timepoint, rep(tp, nrow(phenotypes)))
}
timepoint <- as.numeric(gsub("d", "", timepoint))
topMarker <- as.factor(unlist(genotypes["UNC5048297", individual]))


### Map everything using an across timepoint linear mixed model, using subfamily as a random effect and individual as grouping factor
cl <- makeCluster(getOption("cl.cores", 6))
clusterExport(cl, "phenotype") # export subfamily to the nodes
clusterExport(cl, "individual") # export subfamily to the nodes
clusterExport(cl, "subfamily") # export subfamily to the nodes
clusterExport(cl, "littersize") # export littersize to the nodes
clusterExport(cl, "litternumber") # export litternumber to the nodes
clusterExport(cl, "season") # export season to the nodes
clusterExport(cl, "topMarker") # export season to the nodes
clusterExport(cl, "timepoint") # export toAdd to the nodes
clusterEvalQ(cl, library(lme4)) # load lme4 per node
LMMglobal <- parApply(cl, genotypes[, individual], 1, function(marker){
  noG <- which(is.na(marker))
  if(length(noG) > 0){
    phenotype <- phenotype[-noG]
    individual <- individual[-noG]
    subfamily <- subfamily[-noG]
    littersize <- littersize[-noG]
    litternumber <- litternumber[-noG]
    season <- season[-noG]
    timepoint <- timepoint[-noG]
    topMarker <- topMarker[-noG]
    marker <- marker[-noG]
  }
  subfamily <- as.factor(unlist(subfamily))
  littersize <- as.factor(unlist(littersize))
  litternumber <- as.factor(unlist(litternumber))
  marker <- as.factor(unlist(marker))

  model.fullC   <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) + topMarker*timepoint + marker + marker:timepoint, REML=FALSE)
  model.markerC <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) + topMarker*timepoint + marker, REML=FALSE)
  model.intC    <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) + topMarker*timepoint + marker:timepoint, REML=FALSE)
  model.nullC   <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) + topMarker*timepoint , REML=FALSE)

  model.full    <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) + marker + marker:timepoint, REML=FALSE)
  model.marker  <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) + marker, REML=FALSE)
  model.int     <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) + marker:timepoint, REML=FALSE)
  model.null    <- lmer(phenotype ~ timepoint + littersize + litternumber + season + (1|subfamily/individual) , REML=FALSE)

  
  pvals <- rep(NA, 8)
  names(pvals) <- c("F_M_C", "F_I_C", "M_N_C", "F_N_C", "F_M", "F_I", "M_N", "F_N")

  pvals["F_M_C"] <- anova(model.fullC, model.markerC)["model.full", "Pr(>Chisq)"]
  pvals["F_I_C"] <- anova(model.fullC, model.intC)["model.full", "Pr(>Chisq)"]
  pvals["M_N_C"] <- anova(model.markerC, model.nullC)["model.marker", "Pr(>Chisq)"]
  pvals["F_N_C"] <- anova(model.fullC, model.nullC)["model.full", "Pr(>Chisq)"]
  
  pvals["F_M"] <- anova(model.full, model.marker)["model.full", "Pr(>Chisq)"]
  pvals["F_I"] <- anova(model.full, model.int)["model.full", "Pr(>Chisq)"]
  pvals["M_N"] <- anova(model.marker, model.null)["model.marker", "Pr(>Chisq)"]
  pvals["F_N"] <- anova(model.full, model.null)["model.full", "Pr(>Chisq)"]

  return(pvals)
})
stopCluster(cl)

threshold <- -log10(0.01/nrow(genotypes))

plot(-log10(LMMglobal["F_M", ]), t ='l', col=1)
points(-log10(LMMglobal["F_M_C", ]), t ='l', lwd=2, col='red')
points(20+nA / 10, col='red', pch=19, cex=0.7)
points(20+nH / 10, col='green', pch=19, cex=0.7)
points(20+nB / 10, col='blue', pch=19, cex=0.7)
abline(h = -log10(0.05 / nrow(genotypes)))

threshold <- -log10(0.01/nrow(genotypes))
aboveThreshold <- map[names(which(-log10(LMMglobal["F_M_C", ]) > threshold)),]

regions <- NULL
regionStart <- c(aboveThreshold[1, "Chr"], aboveThreshold[1, "Mb_NCBI38"])
posNow <- c(aboveThreshold[1, "Chr"], aboveThreshold[1, "Mb_NCBI38"])
regionEnd <- c(aboveThreshold[1, "Chr"], aboveThreshold[1, "Mb_NCBI38"])
nInR <- 1
x <- 2
while(x < nrow(aboveThreshold)){
  if(posNow[1] == aboveThreshold[x, "Chr"] && as.numeric(aboveThreshold[x, "Mb_NCBI38"]) < (as.numeric(posNow[2]) + 2500000)){
    posNow <- c(aboveThreshold[x, "Chr"], aboveThreshold[x, "Mb_NCBI38"])
    regionEnd <- c(aboveThreshold[x, "Chr"], aboveThreshold[x, "Mb_NCBI38"])
    nInR <- nInR + 1
    #cat(x, " inR ", posNow[1], posNow[2], aboveThreshold[x, "Chr"], aboveThreshold[x, "Mb_NCBI38"], "\n")
  }else{
    regions <- rbind(regions, c(regionStart, regionEnd, nInR))
    regionStart <- c(aboveThreshold[x, "Chr"], aboveThreshold[x, "Mb_NCBI38"])
    posNow <- c(aboveThreshold[x, "Chr"], aboveThreshold[x, "Mb_NCBI38"])
    regionEnd <- c(aboveThreshold[x, "Chr"], aboveThreshold[x, "Mb_NCBI38"])
    nInR <- 1
    #cat(x, " newR ", posNow[1], posNow[2], aboveThreshold[x, "Chr"], aboveThreshold[x, "Mb_NCBI38"], "\n")
  }
  x <- x +1
}
regions <- rbind(regions, c(regionStart, regionEnd, nInR))
regions <- regions[which(as.numeric(regions[,5]) > 10),]


onChr <- rownames(map)[which(map[,"Chr"] == 3)]
onChr <- onChr[which(onChr %in% rownames(genotypes))]

plot(c(min(map[onChr, "sumPos"]), max(map[onChr, "sumPos"])), c(0, 40), t = 'n', ylab="LOD score", xlab="Chromosome")
for (tp in timepoints) {
  points(map[onChr, "sumPos"], -log10(results.lmm[[tp]]["Full",onChr]), t ='l', col = tpcolors[tp])
}
points(map[onChr, "sumPos"], -log10(LMMglobal["F_M", onChr]) / 8, t ='l', lwd=2, col='blue')
points(map[onChr, "sumPos"], -log10(LMMglobal["F_M_C", onChr]) / 8, t ='l', lwd=2, col='red')
axis(1, at = (chrs.starts + (chrs.lengths / 2))[1:length(chrs)], chrs)
abline(h = -log10(0.05/nrow(genotypes)), lty=2, col="orange")
abline(h = -log10(0.01/nrow(genotypes)), lty=2, col="green")

library(biomaRt)
bio.mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")

kegg <- read.table("mmu04910_kegg.txt",sep='\t', row.names=2, header = FALSE, fill = TRUE)

res.biomart <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "external_gene_name", "mgi_id", "mgi_symbol", "mgi_description"), 
                        filters = c("external_gene_name"), values = rownames(kegg), mart = bio.mart)

res.biomart <- res.biomart[which(res.biomart[,"chromosome_name"] %in% names(chrs)),]
rownames(res.biomart) <- res.biomart[,"external_gene_name"]

kegg <- kegg[which(rownames(kegg) %in% rownames(res.biomart)),]
kegg <- cbind(kegg, res.biomart[rownames(kegg),])
colnames(kegg)[1:4] <- c("kegg_id", "kegg_description", "ko_id", "ec_id")

resultsFULL <- vector("list", nrow(kegg))
names(resultsFULL) <- rownames(kegg)
results <- NULL
for(x in 1:nrow(kegg)){
  moffset <- 0
  region <- c()
  while(length(region) < 2){
    region <- rownames(map[map[,"Chr"] == kegg[x, "chromosome_name"] & 
                      as.numeric(map[,"Mb_NCBI38"]) > (as.numeric(kegg[x, "start_position"]) - moffset) & 
                      as.numeric(map[,"Mb_NCBI38"]) < (as.numeric(kegg[x, "end_position"]) + moffset),])
    moffset <- moffset + 50000
  }
  regiongts <- cbind(nA = nA[region], nH = nH[region], nB = nB[region])
  parentgts <- parentgeno[region,]
  
  resultsFULL[[x]] <- list(moffset = moffset, size = nrow(regiongts), regiongts = regiongts, parentgts = parentgts)
  toSwap <- which(parentgts[,1] == "B")
  if(length(toSwap) > 0){
    temp <- regiongts[toSwap, "nB"]
    regiongts[toSwap, "nB"] <- regiongts[toSwap, "nA"]
    regiongts[toSwap, "nA"] <- temp
  }
  colnames(regiongts) <- c("B6N", "H", "BFMI")
  rGT <- round(apply(regiongts,2, sum, na.rm=TRUE)/nrow(regiongts),2)
  status = "S"
  if(rGT[1] == 0) status = "BFMI"
  if(rGT[3] == 0) status = "B6N"
  results <- rbind(results, c(moffset, nrow(regiongts), rGT, status))
}
rownames(results) <- rownames(kegg)
data.frame(results)



## Chromosome 3 QTL region1:
region1 <- c("UNC4985064", "JAX00105915", "UNC4986011")
map[region1,]
parentgeno[region1,]
## Zmat3, and Pik3ca (growth, Involved in signaling via insulin-receptor substrate (IRS) proteins), (just outside: Kcnmb3 KEGG pathway - Insuline secretion)

## QTL region2, missing homozygous group
region2 <- c("UNC5261852", "UNC5262211","UNC030661179","UNC5263981","JAX00107401")
map[region2,]
parentgeno[region2,]
cbind(nA = nA[region2], nH = nH[region2], nB = nB[region2])
## Foxo1, Major regulator of insuling, -> Never B6N|B6N

#AKT1
akt1mb <- rownames(map[map[,"Chr"] == 12 & as.numeric(map[,"Mb_NCBI38"]) > 111651821 & as.numeric(map[,"Mb_NCBI38"]) < 113676276,])
cbind(nA = nA[akt1mb], nH = nH[akt1mb], nB = nB[akt1mb])

#Ywhaq 12:21390071-21417637:-1 -> Never B6N|B6N
Ywhaqmb <- rownames(map[map[,"Chr"] == 12 & as.numeric(map[,"Mb_NCBI38"]) > 20390071 & as.numeric(map[,"Mb_NCBI38"]) < 22417637,])
cbind(nA = nA[Ywhaqmb], nH = nH[Ywhaqmb], nB = nB[Ywhaqmb])
parentgeno[Ywhaqmb,]

#Ywhae 11:75732869-75765845:1
Ywhaemb <- rownames(map[map[,"Chr"] == 11 & as.numeric(map[,"Mb_NCBI38"]) > 75732869 & as.numeric(map[,"Mb_NCBI38"]) < 75765845,])
cbind(nA = nA[Ywhaemb], nH = nH[Ywhaemb], nB = nB[Ywhaemb])
parentgeno[Ywhaemb,]

#Ywhah 5:33018816-33027966:1
Ywhahmb <- rownames(map[map[,"Chr"] == 5 & as.numeric(map[,"Mb_NCBI38"]) > 32018816 & as.numeric(map[,"Mb_NCBI38"]) < 34027966,])
cbind(nA = nA[Ywhahmb], nH = nH[Ywhahmb], nB = nB[Ywhahmb])
parentgeno[Ywhahmb,]

#Sfn 4:133600556-133602168:-1
Sfnmb <- rownames(map[map[,"Chr"] == 4 & as.numeric(map[,"Mb_NCBI38"]) > 132600556 & as.numeric(map[,"Mb_NCBI38"]) < 134602168,])
cbind(nA = nA[Sfnmb], nH = nH[Sfnmb], nB = nB[Sfnmb])
parentgeno[Sfnmb,]

#Ywhaz 15:36770770-36796929:-1
Ywhazmb <- rownames(map[map[,"Chr"] == 15 & as.numeric(map[,"Mb_NCBI38"]) > 36770770 & as.numeric(map[,"Mb_NCBI38"]) < 36796929,])
cbind(nA = nA[Ywhazmb], nH = nH[Ywhazmb], nB = nB[Ywhazmb])
parentgeno[Ywhazmb,]

#Ywhab 2:163994960-164018588:1 -> Never BFMI|BFMI
Ywhabmb <- rownames(map[map[,"Chr"] == 2 & as.numeric(map[,"Mb_NCBI38"]) > 162994960 & as.numeric(map[,"Mb_NCBI38"]) < 164994960,])
cbind(nA = nA[Ywhabmb], nH = nH[Ywhabmb], nB = nB[Ywhabmb])
parentgeno[Ywhabmb,]


## QTL region3, missing homozygous group
region3 <- c("UNC5475952","JAX00526202","UNC5482456","UNC5486706","UNC5492395","UNC5496983","UNC5514046","backupUNC030098041","UNC5532424","UNC030101822", "UNC5560335")
map[region3,]
## Bche (Subcutaneous fat pad expression), Sis, Slitrk3, Serpini1 and Serpini 2
## Butyrylcholinesterase Deficiency Promotes Adipose Tissue Growth and Hepatic Lipid Accumulation in Male Mice on High-Fat Diet.
