getSeason <- function(mmonths) {
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5] <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8] <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11] <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2] <- "Winter"
  return(ret)
}

setwd("D:/Edrive/Cow/PaulaGWAS")
load("inputGWAS_1800DSNcows.Rdata")

father <- pdata[,"father"]
farm <- pdata[,"farm"]
birthyear <- substr(pdata[,"D_GEB"], 0, 4)
kalfyear <- substr(pdata[,"D_KALB_1"], 0, 4)
kalfseason <- getSeason(as.numeric(substr(pdata[,"D_KALB_1"], 5, 6)))
mkg <- pdata[,"LA_MKG_1"]

mdata <- data.frame(mkg, farm, father, birthyear, kalfyear, kalfseason)
rownames(mdata) <- rownames(pdata)
mdata <- mdata[-which(apply(mdata,1, function(x){any(is.na(x))})),]

# Cows need to have a fathers with more than 20 offspring
mdata <- mdata[which(pdata[,"father"] %in% names(which(table(mdata[,"father"]) > 20))),]

# Cows need to be from a year with more than 20 observations
mdata <- mdata[which(substr(mdata[,"birthyear"], 0, 4) %in% names(which(table(substr(mdata[,"birthyear"], 0, 4)) > 20))),]

# Cows need to be from a farm with more than 20 observations
mdata <- mdata[which(mdata[,"farm"] %in% names(which(table(mdata[,"farm"]) > 20))),]

dim(mdata)

genotypes <- genotypes[, rownames(mdata)]

gtsnew <- t(apply(genotypes,1, function(x){
  tbl <- table(x)
  if(any(tbl < 30)){
    x[which(x == names(which(tbl < 30)))] <- NA
  }
  return(x)
}))

library(parallel)

cl <- makeCluster(6)
pvalM <- t(parApply(cl, gtsnew, 1, function(m, mdata){
  add <- m
  dd <- as.numeric(m == 0)
  ret <- rep(NA, 7)
  res <- anova(lm(mkg ~ as.factor(farm) + as.factor(father) + as.factor(birthyear) + as.factor(kalfyear) + as.factor(kalfseason) + add + dd, data = mdata))
  ret[1:(length(res[[5]])-1)] <- res[[5]][1:(length(res[[5]])-1)]
  return(ret)
}, mdata = mdata))
stopCluster(cl)

observed.p.add <- pvalM[,6]
expected.p.add <- (rank(observed.p.add, ties.method="first")+0.5) / (length(observed.p.add) + 1)
expected.l.add <- -log10(sort(expected.p.add, na.last = TRUE))
observed.l.add <- -log10(sort(observed.p.add, na.last = TRUE))

round(median(qchisq(1.0 - observed.p.add, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)