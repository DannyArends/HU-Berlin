MDCfamily <- "7984"

setwd("D:/Edrive/Mouse/MDC/Jan2021")
geno <- read.csv(paste0("geno_", MDCfamily, ".txt"), sep = "\t", na.strings = c("", "NA", "-", "X"))

pheno <- read.csv(paste0("pheno_", MDCfamily, ".txt"), sep = "\t", na.strings = c("", "NA", "-", "X", "#DIV/0!"))
pheno2 <- read.csv(paste0("pheno_", MDCfamily, "_mri.txt"), sep = "\t", na.strings = c("", "NA", "-", "X", "#DIV/0!"))
pheno2 <- pheno2[-which(is.na(pheno2[, "ID.Nr"])),]

colnames(geno) <- c("MouseID", "GT", "Uncertain")

dim(pheno)
dim(geno)

mother <- NA
father <- NA
wdate <- NA
WG <- NA
GEN <- NA
for(x in 1:nrow(pheno)){
  if(!is.na(pheno[x, "Vater"])){
    mother <- pheno[x, "Mutter"]
    father <- pheno[x, "Vater"]
    wdate <- pheno[x, "W.dat"]
    WG <- pheno[x, "WG"]
    GEN <- pheno[x, "Gen."]
  }else{
    pheno[x, "Mutter"] <- mother
    pheno[x, "Vater"] <- father
    pheno[x, "W.dat"] <- wdate
    pheno[x, "WG"] <- WG
    pheno[x, "Gen."] <- GEN
  }
}

sex <- as.character(pheno[,"m"])
sex[which(is.na(sex))] <- "f"
pheno <- cbind(pheno, sex = sex)

pheno <- pheno[, c("ID.Nr", "sex", "Gen.", "Mutter", "Vater", "W.dat", "WG", "X21", "X28", "X35", "X42", "X49", "X56", "X63", "X70")]

pheno <- cbind(pheno, GenGood = NA)
lookupTable <- c()
for(x in 1:nrow(pheno)){
  mother <- gsub("MDC-", "", as.character(pheno[x, "Mutter"]))
  father <- gsub("MDC-", "", as.character(pheno[x, "Vater"]))
  child <- gsub("MDC-", "", as.character(pheno[x, "ID.Nr"]))
  
  mGen <- NA
  fGen <- NA
  cGen <- NA

  if(substr(mother, 0, 3) == "100") mGen <- 0
  if(substr(mother, 0, 3) == "101") mGen <- 0
  if(substr(mother, 0, 3) == "TCF") mGen <- 0
  if(substr(mother, 0, 3) == "860") mGen <- 0
  if(substr(mother, 0, 3) == "816") mGen <- 0
  
  if(substr(father, 0, 3) == "100") fGen <- 0
  if(substr(father, 0, 3) == "101") fGen <- 0
  if(substr(father, 0, 3) == "TCF") fGen <- 0
  if(substr(father, 0, 3) == "860") fGen <- 0
  if(substr(father, 0, 3) == "816") fGen <- 0
  
  if(is.na(mGen)){
    #Figure out what generation the mother is
    idx <- which(lookupTable[,1] == mother)
    if(length(idx == 1)){
      mGen <- as.numeric(lookupTable[idx, 2])
    }else{
      stop(paste0("mother ",mother," not in lookup table"))
    }
  }
  if(is.na(fGen)){
    #Figure out what generation the father is
    idx <- which(lookupTable[,1] == father)
    if(length(idx == 1)){
      fGen <- as.numeric(lookupTable[idx, 2])
    }else{
      stop(paste0("father ",father," not in lookup table"))
    }
  }
  if(!is.na(mGen) && !is.na(fGen) && mGen == 0 && fGen == 0){
    cGen <- 1
    lookupTable <- rbind(lookupTable, c(child, cGen))
  }
  if(!is.na(mGen) && !is.na(fGen) && max(mGen, fGen) > 0){
    cGen <- max(mGen, fGen) + 1
    lookupTable <- rbind(lookupTable, c(child, cGen))
  }
  pheno[x, "GenGood"] <- cGen
}

# set the name of the animal as the rowname !!!
rownames(pheno) <- gsub("MDC-", "", pheno[, "ID.Nr"])
rownames(pheno2) <- gsub("MDC-", "", pheno2[, "ID.Nr"])
pheno2 <- pheno2[which(rownames(pheno2) %in% rownames(pheno)),]

pheno <- cbind(pheno, pheno2[rownames(pheno),c("Fat.56", "Fat.70", "Fat.Lean56", "Fat.Lean70", "FAT56", "FAT70", "LEAN56", "LEAN70")])

dupMouseID <- geno[duplicated(geno[, "MouseID"]),"MouseID"]

for(m in dupMouseID){
  dupIDs <- which(geno[, "MouseID"] == m)
  tbl <- table(geno[dupIDs, "GT"])
  if(length(tbl == 1)){
    # Duplicates have same genotypes
    if(any(!is.na(geno[dupIDs, "Uncertain"]))){
      cat("Uncertainty for mouse: ",m,"\n")
      notUncertain <- which(is.na(geno[dupIDs, "Uncertain"]))
      geno <- geno[-dupIDs[-notUncertain],]
    }else{
      geno <- geno[-dupIDs[-1],]
    }
  }else{
    cat("Multiple GTS for mouse: ",m,"\n")
    geno <- geno[-dupIDs,]
  }
}

rownames(geno) <- paste0(MDCfamily,"-", geno[, "MouseID"])

mdata <- cbind(pheno, geno[rownames(pheno),])

mdata <- mdata[-which(is.na(mdata[, "GT"])),]
if(length(which(!is.na(mdata[, "Uncertain"]))) > 0) mdata <- mdata[-which(!is.na(mdata[, "Uncertain"])),]
mdata <- mdata[, -which(colnames(mdata) %in% c("Uncertain", "Ignore"))]

#hasNA <- which(apply(apply(mdata, 1, is.na),2,sum) > 0)
#if(length(hasNA) > 0){
#  mdata <- mdata[-hasNA,]
#}

# Correct for small genotype groups
smallGTS <- which(mdata[, "GT"] %in% names(which(table(mdata[, "GT"]) < 10)))
if(length(smallGTS) > 0){ mdata <- mdata[-smallGTS,] }

correctPhe <- function(mdata, column = "X21", name = "C21"){
  cresiduals <- residuals(lm(as.numeric(mdata[, column]) ~ mdata[, "sex"] + as.factor(mdata[, "WG"])))
  cres <- rep(NA, nrow(mdata))
  cres[as.numeric(names(cresiduals))] <- cresiduals
  corrected <- round(mean(as.numeric(mdata[, column]),na.rm=TRUE) + cres,2)
  mdata <- cbind(mdata, corrected)
  colnames(mdata)[ncol(mdata)] <- name
  return(mdata)
}
dim(mdata)
mdata <- correctPhe(mdata, "X21", "C21")
mdata <- correctPhe(mdata, "X21", "C28")
mdata <- correctPhe(mdata, "X35", "C35")
mdata <- correctPhe(mdata, "X42", "C42")
mdata <- correctPhe(mdata, "X49", "C49")
mdata <- correctPhe(mdata, "X56", "C56")
mdata <- correctPhe(mdata, "X63", "C63")
mdata <- correctPhe(mdata, "X70", "C70")
mdata <- correctPhe(mdata, "Fat.Lean56", "FatLeanC56")
mdata <- correctPhe(mdata, "Fat.Lean70", "FatLeanC70")

for(x in c("C21","C35","C42","C49","C56","C63","C70")){
  minV <- mean(mdata[,x],na.rm=TRUE) - (3 * sd(mdata[,x],na.rm=TRUE))
  maxV <- mean(mdata[,x],na.rm=TRUE) + (3 * sd(mdata[,x],na.rm=TRUE))
  mdata[which(mdata[,x] < minV | mdata[,x] > maxV), x] <- NA
}

GTS <- rep(NA, nrow(mdata))
GTS[mdata[, "GT"] == "11"] <- "BFMI/BFMI"
GTS[mdata[, "GT"] == "12"] <- "BFMI/B6N"
GTS[mdata[, "GT"] == "33"] <- "MDC/MDC"
GTS[mdata[, "GT"] == "31"] <- "MDC/BFMI"
GTS[mdata[, "GT"] == "13"] <- "MDC/BFMI"
GTS[mdata[, "GT"] == "32"] <- "MDC/B6N"
GTS[mdata[, "GT"] == "23"] <- "MDC/B6N"
GTS[mdata[, "GT"] == "22"] <- "B6N/B6N"

mdata[, "GT"] <- GTS

#hasNA <- which(apply(apply(mdata, 1, is.na),2,sum) > 0)
#if(length(hasNA) > 0){ mdata <- mdata[-hasNA,] }

phe <- "C70"

dim(mdata)
plot(x = c(0.5, 0.5+length(unique(mdata[,"GT"]))), y = c(0,0.70), t = 'n', xaxs="i", yaxs="i", xaxt='n', las=2, xlab="GT", ylab=paste0(phe))
mybplot <- boxplot(mdata[, phe] ~ mdata[,"GT"], main = paste0("Family = ", MDCfamily, " (", phe, ")"), add = TRUE, yaxt='n', xaxt='n')
nums <- c()
for(gt in mybplot$names){
  n <- length(which(!is.na(mdata[which(mdata[,"GT"] == gt), phe])))
  nums <- rbind(nums, c(gt, n))
}
xn <- apply(nums,1, paste0,collapse="\n n=")
axis(1, at = 1:length(xn), xn)

mdata[which(mdata[,"GT"] == "MDC/BFMI"), phe]
mdata[which(mdata[,"GT"] == "MDC/B6N"), phe]
mdata[which(mdata[,"GT"] == "BFMI/BFMI"), phe]

for(gt in unique(mdata[,"GT"])){
  cat(gt, mean(mdata[which(mdata[,"GT"] == gt), phe],na.rm=TRUE), "+/-", sd(mdata[which(mdata[,"GT"] == gt), phe],na.rm=TRUE), "\n")
}


p12 <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/B6N", phe], mdata[mdata[,"GT"] == "BFMI/BFMI", phe])$p.value},error = function(e) {return(NA)})
p13 <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/B6N", phe], mdata[mdata[,"GT"] == "MDC/B6N", phe])$p.value},error = function(e) {return(NA)})
p14 <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/B6N", phe], mdata[mdata[,"GT"] == "MDC/BFMI", phe])$p.value},error = function(e) {return(NA)})
p15 <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/B6N", phe], mdata[mdata[,"GT"] == "MDC/MDC", phe])$p.value},error = function(e) {return(NA)})

p23 <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/BFMI", phe], mdata[mdata[,"GT"] == "MDC/B6N", phe])$p.value},error = function(e) {return(NA)})
p24 <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/BFMI", phe], mdata[mdata[,"GT"] == "MDC/BFMI", phe])$p.value},error = function(e) {return(NA)})
p25 <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/BFMI", phe], mdata[mdata[,"GT"] == "MDC/MDC", phe])$p.value},error = function(e) {return(NA)})

p34 <- tryCatch({t.test(mdata[mdata[,"GT"] == "MDC/B6N", phe], mdata[mdata[,"GT"] == "MDC/BFMI", phe])$p.value},error = function(e) {return(NA)})
p35 <- tryCatch({t.test(mdata[mdata[,"GT"] == "MDC/B6N", phe], mdata[mdata[,"GT"] == "MDC/MDC", phe])$p.value},error = function(e) {return(NA)})

p45 <- tryCatch({t.test(mdata[mdata[,"GT"] == "MDC/BFMI", phe], mdata[mdata[,"GT"] == "MDC/MDC", phe])$p.value},error = function(e) {return(NA)})

asN <- function(x){
  return(formatC(x, format = "e", digits = 1))
}

asNN <- function(x){
  format(round(x, 2), nsmall = 2)
}

colz <- rep("black",10)
colz[3:4] <- "red"
colz[8:9] <- "darkgreen"

options(digits = 2)
legend("topright", 
  c(paste0("BFMI/B6N v BFMI/BFMI (p < 0.05)"),#", asN(min(p12 * 10,1)), ")"), 
    paste0("BFMI/B6N v MDC/B6N (N.S.)"),#, asNN(min(p13 * 10,1)), ")"),
    paste0("BFMI/B6N v MDC/BFMI (p = 0.58)"),#, asNN(min(p14 * 10,1)), ")"),
    paste0("BFMI/B6N v MDC/MDC (p = 0.77)"),#, asNN(min(p15 * 10,1)), ")"),
    paste0("BFMI/BFMI v MDC/B6N (p < 0.05)"),#, asN(min(p23 * 10,1)), ")"),
    paste0("BFMI/BFMI v MDC/BFMI (p < 0.05)"),#, asN(min(p24 * 10,1)), ")"),
    paste0("BFMI/BFMI v MDC/MDC (p < 0.05)"),#, asN(min(p25 * 10,1)), ")"),
    paste0("MDC/B6N v MDC/BFMI (p < 0.05)"),#, asN(min(p34 *10,1)), ")"),
    paste0("MDC/B6N v MDC/MDC (p < 0.05)"),#, asN(min(p35*10,1)), ")"),
    paste0("MDC/BFMI v MDC/MDC (N.S.)"))#, asNN(min(p45*10,1)), ")"))
    , text.col = colz)



BvH <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/BFMI", phe], mdata[mdata[,"GT"] == "MDC/BFMI", phe])$p.value},error = function(e) {return(NA)})
BvM <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/BFMI", phe], mdata[mdata[,"GT"] == "MDC/MDC", phe])$p.value},error = function(e) {return(NA)})
HvM <- tryCatch({t.test(mdata[mdata[,"GT"] == "MDC/BFMI", phe], mdata[mdata[,"GT"] == "MDC/MDC", phe])$p.value},error = function(e) {return(NA)})
BvB <- tryCatch({t.test(mdata[mdata[,"GT"] == "BFMI/BFMI", phe], mdata[mdata[,"GT"] == "B6N/B6N", phe])$p.value},error = function(e) {return(NA)})
HvB <- tryCatch({t.test(mdata[mdata[,"GT"] == "MDC/BFMI", phe], mdata[mdata[,"GT"] == "B6N/B6N", phe])$p.value},error = function(e) {return(NA)})
MvB <- tryCatch({t.test(mdata[mdata[,"GT"] == "MDC/MDC", phe], mdata[mdata[,"GT"] == "B6N/B6N", phe])$p.value},error = function(e) {return(NA)})


legend("topright", 
  c(paste0("BFMI/BFMI v MDC/BFMI (p = ", round(BvH,4), ")"), 
    paste0("BFMI/BFMI v MDC/MDC (p = ", round(BvM,4), ")"),
    paste0("MDC/BFMI v MDC/MDC (p = ", round(HvM,4), ")"),
    paste0("BFMI/BFMI v B6N/B6N (p = ", round(BvB,4), ")"),
    paste0("MDC/BFMI v B6N/B6N (p = ", round(HvB,4), ")"),
    paste0("MDC/MDC v B6N/B6N (p = ", round(MvB,4), ")")))


#complementation <- mdata[which(mdata[,"GT"] == "MDC/BFMI"), phe]
#mrange <- seq(5, 80, 5)
#nps <- c()
#for(x in mrange){
#  np <- 0
#  for(p in 1:1000){
#    null <- rnorm(x, 0.096, 0.036)
#    if(t.test(complementation, null)$p.value < 0.005) np <- np+1
#  }
#  nps <- c(nps, np / 1000)
#}
#plot(nps, xaxt='n', ylab ="Power", xlab="N samples in negative control", pch=19, las=2)
#axis(1, at = 1:length(mrange), mrange)
#abline(h = 0.80, col="red")
#abline(h = 0.95, col="orange")
#abline(h = 0.995, col = "green")
#legend("bottomright", c("80% power (original estimate)", "95% (no multiple test adjustment)", "99% (multiple test adjusted)"), lwd=1, col=c("red", "orange", "green"))


# too many generations only plot the last 8
# 
op <- par(mfrow=c(2,3))  
for(g in unique(mdata[, "Gen."])[(length(unique(mdata[, "Gen."]))-5): length(unique(mdata[, "Gen."]))]){
  cdata <- mdata[which(mdata[, "Gen."] == g),]
  plot(x = c(0.5, 0.5+length(unique(cdata[,"GT"]))), y = c(0,0.50), t = 'n', 
       xaxs="i", yaxs="i", xaxt='n', las=2, xlab="GT", ylab=phe, main = paste0("Family = ", MDCfamily, ", Gen = ", g))
  mybplot <- boxplot(cdata[, phe] ~ cdata[,"GT"], add = TRUE, yaxt='n', xaxt='n')
  nums <- c()
  for(gt in mybplot$names){
    n <- length(which(!is.na(cdata[which(cdata[,"GT"] == gt), phe])))
    nums <- rbind(nums, c(gt, n))
  }
  xn <- apply(nums,1, paste0,collapse="\n n=")
  axis(1, at = 1:length(xn), xn)
}


