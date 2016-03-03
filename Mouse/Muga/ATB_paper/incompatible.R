incompatible <- function(marker1, marker2, homozygous = FALSE){
  nodata <- c(which(is.na(marker1)), which(is.na(marker2)))
  if(length(nodata) > 0){
    marker1 <- marker1[-nodata]
    marker2 <- marker1[-nodata]
  }
  mR1 <- table(as.character(marker1))
  mR2 <- table(as.character(marker2))

  combinations <- NULL
  for(a1 in names(mR1)){
    for(a2 in names(mR2)){
      combinations <- c(combinations, paste0(a1, a2))
    }
  }

  observed <- rep(0, length(combinations)); names(observed) <- combinations
  for(i in paste0(marker1,marker2)){
    observed[i] <- observed[i] + 1
  }
  expected <- rep(0, length(combinations)); names(expected) <- combinations
  for(a1 in names(mR1)){
    for(a2 in names(mR2)){
      expected[paste0(a1, a2)] <- (mR1[a1] / sum(mR1)) * (mR2[a2] / sum(mR2))
    }
  }
  expected <- round(sum(observed) * expected, 0)
  
  chiSQ <- 0
  for(combi in combinations) {
    if(homozygous && grepl("H", combi)){
      # Do nothing for heterozygous individuals
    }else{
      if(expected[combi] != 0){
        n <- (as.numeric(observed[combi] - expected[combi]) ^ 2) / expected[combi]
        chiSQ <- chiSQ + n
      }
    }
  }

  list(obs = observed, exp = expected, chisq = chiSQ)
}

doLD <- function(marker1, marker2){
  nodata <- c(which(is.na(marker1)), which(is.na(marker2)))
  if(length(nodata) > 0){
    marker1 <- marker1[-nodata]
    marker2 <- marker1[-nodata]
  }
  m1o <-  rep("", length(marker1))
  m1o[marker1 == "A"] <- "A/A"
  m1o[marker1 == "H"] <- "A/T"
  m1o[marker1 == "B"] <- "T/T"

  m2o <-  rep("", length(marker2))
  m2o[marker2 == "A"] <- "A/A"
  m2o[marker2 == "H"] <- "A/T"
  m2o[marker2 == "B"] <- "T/T"

  unlist(LD(genotype(m1o),genotype(m2o)))$"X^2"
}

MAF <- function(marker){
  nodata <- which(is.na(marker))
  if(length(nodata) > 0) marker <- marker[-nodata]

  mR1 <- table(as.character(marker))
  minAllele <- names(which.min(mR1))
  return(as.numeric((mR1[minAllele] *2 + mR1["H"]) / (2*sum(mR1))))
}

generate <- function(nind = 100, nrecom = 1) {
  m1 <- sample(c("A","H","B"), nind, replace=TRUE)
  m2 <- m1
  for(x in 1:nrecom){
    l <- sample(1:nind, 1)
    if(m2[l] == "A") m2[l] <- "H"
    if(m2[l] == "H"){
      if(runif(1) < 0.5){ 
        m2[l] <- "A"
      }else{
        m2[l] <- "B"
      }
    }
    if(m2[l] == "B") m2[l] <- "H"
  }
  list(m1,m2)
}

library(genetics)

resultsMine <- NULL
resultsLD <- NULL
resultsMAFS <- NULL
for(r in 1:100){
  scoresMine <- NULL
  scoresLD <- NULL
  scoresMAFs <- NULL
  for(x in 1:10){
    markers <- generate(50, r)
    scoresMine <- c(scoresMine, incompatible(markers[[1]],markers[[2]])$chisq)
    scoresLD <- c(scoresLD, doLD(markers[[1]],markers[[2]]))
    scoresMAFs <- c(scoresMAFs, MAF(markers[[1]]) - MAF(markers[[2]]))
  }
  resultsMine <- rbind(resultsMine, scoresMine)
  resultsLD <- rbind(resultsLD, scoresLD)
  resultsMAFS <- rbind(resultsMAFS, scoresMAFs)
  cat(r, "\n")
}

boxplot(t(resultsMine))
boxplot(t(resultsLD), col="red", add=TRUE)

plot(resultsMine[,1], resultsLD[,1])
cor(resultsMine[,1], resultsLD[,1])