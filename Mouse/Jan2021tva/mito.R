# Analysis of mito line 1 versus pure BFMI
setwd("D:/Edrive/Mouse/BFMI860-12")
mdpure <- read.table("purebred.txt",sep="\t", skip=1, header=TRUE,na.strings=c("NA", "WFV", ""))
mdmito <- read.table("mito.txt",sep="\t", skip=2, header=TRUE)

for(x in 1:nrow(mdmito)){
  if(is.na(mdmito[x, "WG"])){ mdmito[x, "WG"] <- WG;
  }else{ mdmito[x, "WG"] -> WG; }
}

for(x in 1:nrow(mdpure)){
  if(is.na(mdpure[x, "WG21"])){ mdpure[x, "WG21"] <- WG;
  }else{ mdpure[x, "WG21"] -> WG; }
}

iiP <- which(mdpure[, "WG21"] >= 7 & mdpure[, "WG21"] <= 9 & mdpure[, "Sex"] == "m")
iiM <- which(mdmito[, "WG"] >= 7 & mdmito[, "WG"] <= 9 & mdmito[, "Sex"] == "m")

means.male <- c()
sds.male <- c()
p.male <- c()
for(x in c("X21", "X28", "X35", "X42", "X49", "X56")){
  mdP <- na.omit(as.numeric(mdpure[iiP,x]))
  mdM <- na.omit(as.numeric(mdmito[iiM,x]))
  cat("pure=", mean(mdP), "n=", length(mdP), "mito=", mean(mdM), "n=", length(mdM), "\n")
  means.male  <- rbind(means.male, c(mean(mdP), mean(mdM)))
  sds.male  <- rbind(sds.male, c(sd(mdP), sd(mdM)))
  p.male  <- c(p.male, t.test(mdP, mdM)$p.value)
}
colnames(means.male) <- colnames(sds.male) <- c("BFMI", "MITO")

iiP <- which(mdpure[, "WG21"] >= 7 & mdpure[, "WG21"] <= 9 & mdpure[, "X"] == "f")
iiM <- which(mdmito[, "WG"] >= 7 & mdmito[, "WG"] <= 9 & mdmito[, "X"] == "f")

means.female <- c()
sds.female <- c()
p.female <- c()
for(x in c("X21", "X28", "X35", "X42", "X49", "X56")){
  mdP <- na.omit(as.numeric(mdpure[iiP,x]))
  mdM <- na.omit(as.numeric(mdmito[iiM,x]))
  cat("pure=", mean(mdP), "n=", length(mdP), "mito=", mean(mdM), "n=", length(mdM), "\n")
  means.female  <- rbind(means.female, c(mean(mdP), mean(mdM)))
  sds.female  <- rbind(sds.female, c(sd(mdP), sd(mdM)))
  p.female  <- c(p.female, t.test(mdP, mdM)$p.value)
}
colnames(means.female) <- c("BFMI", "MITO")
colnames(means.female) <- colnames(sds.female) <- c("BFMI", "MITO")

tps <- c(21, 28, 35, 42, 49, 56)

plot(c(20, 60), c(0, 40), t = 'n', xlab="Age (days)", ylab="Bodyweight (g)", las=2)
for(x in 1:nrow(means.male)){
  points(c(tps[x]-1, tps[x]+1), means.male[x,], pch = c(16, 21), col = "blue")
  segments(tps[x]-1, means.male[x,1] - sds.male[x,1], tps[x]-1, means.male[x,1] + sds.male[x,1], col = "blue")
  segments(tps[x]+1, means.male[x,2] - sds.male[x,2], tps[x]+1, means.male[x,2] + sds.male[x,2], col = "blue")

  points(c(tps[x]-0.5, tps[x]+0.5), means.female[x,], pch = c(16, 21), col = "pink")
  segments(tps[x]-0.5, means.female[x,1] - sds.female[x,1], tps[x]-0.5, means.female[x,1] + sds.female[x,1], col = "pink")
  segments(tps[x]+0.5, means.female[x,2] - sds.female[x,2], tps[x]+0.5, means.female[x,2] + sds.female[x,2], col = "pink")
}
legend("bottomright", c("Male(BFMI)", "Male(Mito)", "Female(BFMI)", "Female(Mito)"), col=c("blue","blue","pink","pink"), pch=c(16,21,16,21))

cbind(means.male, p.male)
cbind(means.female, p.female)


calculatepower <- function(m1 = 1.41, m2 = 1.81, sd1=0.66, sd2=sd1, alpha = 0.05, from = 10, to = 100, iterations = 100, alternative = c("less", "two.sided"), verbose = TRUE){
  if(from < 3)  stop("Need at least 3 samples for a t.test")
  if(to < from) stop("to needs to be larger then from:", to, "<", from)
  npoints <- 100
  mvals <- NULL
  for(s in from:to){                                              # s = number of samples
    power <- NULL
    for(n in 1:npoints){                                          # n = replicate to create a boxplot, we fix this to doing 10
      pvals <- NULL
      for(x in 1:iterations){                                     # x = number of times to perform the test (more = better)
        ttest <- t.test(rnorm(s, m=m1, sd=sd1), rnorm(s, m=m2, sd=sd2), alternative = alternative[1])
        pvals <- c(pvals, ttest$p.value)
      }
      power <- c(power, sum(pvals < alpha) / iterations)
    }
    if(verbose) cat("Finished with",s,"\n")
    mvals <- cbind(mvals, power)
  }
  colnames(mvals) <- from:to
  rownames(mvals) <- 1:npoints
  boxplot(mvals, xlab="Number of samples per group", ylab="Power")
  return(invisible(mvals))
}

getpower <- function(powercurve, power = 0.80){
  lowestpower <- rev(apply(powercurve, 2, mean) - apply(powercurve, 2, sd))
  for(x in 1:(length(lowestpower)-1)){
    if(lowestpower[x] > power && lowestpower[x+1] < power){
      cat("You will need in each group:", names(lowestpower)[x], "animals\n")
      return(invisible(names(lowestpower)[x]))
    }
  }
  cat("Unable to find a number of animals, please visually inspect the powercurve\n")
  return(invisible(NULL))
}

powercurveMale <- calculatepower(18.31687, 19.55792, 2.935155, 2.452258, alternative = "less")
powercurveFemale <- calculatepower(17.17358, 18.20457, 2.518310, 2.050462, alternative = "two.sided")
getpower(powercurveMale, 0.80)
getpower(powercurveFemale, 0.80)



boxplot(powercurveMale)