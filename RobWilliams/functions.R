#
# Find recombinations in the most general way possible (by index or location)
# copyright (c) 2016-2020 - Danny Arends and Rob Williams
# last modified Apr, 2016
# first written Apr, 2016
#

recombinations <- function(x, locs = 1:length(x), chr = rep(1, length(x)), fancy = FALSE) {
  if(length(na.omit(x)) < 2) stop("Input x must have at least 2 elements")
  cindex <- 1; cchar <- x[cindex]                 #cchar is the first element
  lindex <- 1
  recombs <- NULL
  while(cindex  <= length(x)) {
    if(is.na(cchar)) {
      cchar <- x[cindex]
    } else { # lchar is now the last non-na character we seen
      if(!is.na(x[cindex])) {
        if(cchar != x[cindex]) {
          location <- (locs[cindex] + locs[lindex]) / 2
          spanschromosomes <- chr[cindex] != chr[lindex]
          if(fancy) location <- paste0(chr[cindex], ":", location)
          if(!spanschromosomes) recombs <- c(recombs, location)
          cchar <- x[cindex]
        }
        lindex <- cindex
      }
    }
    cindex = cindex + 1
  }
  return(recombs)
}

get.locs <- function(x){
  unlist(lapply(lapply(lapply(x, strsplit, ":"), unlist),"[",2))
}

get.chr <- function(x){
  unlist(lapply(lapply(lapply(x, strsplit, ":"), unlist),"[",1))
}

recombinations(c("A","A","A","A","B"))
recombinations(c("A","B","A","A","B"))
recombinations(c("A","B", NA,"A","B"))
recombinations(c("A","B","A","A","B"), c(5, 20, 27, 45, 50) )
recombinations(c("A","B", NA,"A","B"), c(5, 20, 27, 45, 50) )
recombinations(c("A","B", NA,"A","B"), c(5, 20, 27, 45, 50) , c(1,1,2,2,2))
recombinations(c("A","B", NA,"A","B"), c(5, 20, 27, 45, 50) , c(1,1,2,2,2), TRUE)
recombinations(c("A","B", NA,"A","B"), chr = c(1,1,2,2,2), fancy = FALSE)
