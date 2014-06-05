# powerCalculationTtest.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Functions to analyse the required sample size in t.tests

calculatepower <- function(m1 = 1.41, m2 = 1.81, sd1=0.66, sd2=sd1, alpha = 0.05, from = 3, to = 100, iterations = 100, alternative = c("less", "two.sided"), verbose = TRUE){
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

powercurve <- calculatepower()
getpower(powercurve, 0.80)