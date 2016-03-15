tossCoin <- function(){
  if(runif(1) < 0.5) return("H")
  return("T")
}

throwUntill <- function(seq = c("H", "H")){
  res <- NULL
  for(x in 1:length(seq)){
    res <- c(res, tossCoin())
  }
  while(!all(res[(length(res) - length(seq) + 1): length(res)] == seq)){
    res <- c(res, tossCoin())
  }
  return(length(res))
}

throws <- function(n = 2){
  res <- NULL
  for(x in 1:n){
    res <- c(res, tossCoin())
  }
  return(paste0(res,collapse=""))
}

alice <- NULL
bob <- NULL
for(x in 1:10000){
  alice <- c(alice, throwUntill(c("H", "T")))
  bob <- c(bob, throwUntill(c("H", "H")))
}

res <- NULL
for(x in 1:10000){
  res <- c(res, throws(3))
}
sum(table(res)[which(grepl("HT", names(table(res))))]) / sum(table(res))
sum(table(res)[which(grepl("TT", names(table(res))))]) / sum(table(res))