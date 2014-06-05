# \file Answers - Making functions II.R
#
# Copyright (c) 2014-2016, Danny Arends
# Last modified:		June, 2014
# Created:		      May, 2014
# 
# A copy of the GNU General Public License, version 3, is available
# at http://www.r-project.org/Licenses/GPL-3
# 
# Contains: Answers to the Making functions II practical

setwd("~/practicals/Making functions II") # Set working dir

# 0)
factorNR <- function(x){
  fact <- 1
  while(x > 1){
    fact = fact * x
    x = x - 1
  }
  return(fact)
}

factorR <- function(x){
  if(x == 1) return(1)
  return(x * factorR(x - 1))
}

timesNR <- NULL
timesR <- NULL
for(x in 1:1000){
  s <- proc.time()
  factorNR(x)
  nr <- proc.time()
  timesNR <- c(timesNR, (nr-s)[3])
  factorR(x)
  r <- proc.time()
  timesR <- c(timesR, (r-nr)[3])
}

# 1)

fibonacci <- function(x){
  if(x == 0) return(0)
  if(x == 1) return(1)
  return(fibonacci(x-1) + fibonacci(x-2))
}

fibs <- 0
for(x in 1:20){
  fibs <- c(fibs, fibonacci(x))
}

# 2)

Lapply <- function(x, FUN){
  if(class(x) != "list") stop("You need to provide a list to Lapply")
  ret <- x
  for(element in 1:length(x)){
    ret[[element]] <- FUN(x[[element]])
  }
  return(ret)
}

# 3)

F <- function(x){ return(x+1); }
G <- function(x){ return(x/2); }

composition <- function(x, selector, F, G){
  if(selector < 0 || selector > 4) stop("No such selection possible")
  if(selector == 1) return(F(x))
  if(selector == 2) return(G(x))
  if(selector == 3) return(F(G(x)))
  if(selector == 4) return(G(F(x)))
}

# 4)

rawdata <- read.table("statistics.txt")

anova(lm(rawdata[,"Survival"] ~ rawdata[,"Treatment"]))
anova(lm(rawdata[,"Survival"] ~ rawdata[,"Age"]))
anova(lm(rawdata[,"Survival"] ~ rawdata[,"Sex"]))


anova(lm(rawdata[,"Survival"] ~ rawdata[,"Sex"] + rawdata[,"Treatment"]))

# Answers: Treatment, Sex, They shouldn't but they do, Seems so, Yes

# AA1) Euclidean algorithm.

GCD <- function(a, b){
  if(a == b) return(a)
  if(a > b) return(GCD(a-b, b))
  if(b > a) return(GCD(a, b-a))
}

# AA2)
isPrime <- function(x){
  i <- 2
  while(i >= 2 && i <= x / 2){
    if(x %% i == 0) return(FALSE)
    i <- i + 1
  }
  return(TRUE)
}

checkPrimes <- function(){
  primes <- NULL
  for(x in 1:100){
    if(isPrime(x)) primes <- c(primes, x)
  }
  return(primes)
}

# Yes and No, there are partial recursive algorithms, but none that generates all the primes

