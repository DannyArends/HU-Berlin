
toVCFcode <- function(x){
  r <- rep(NA,length(x));  r[x == "A"] <- "0/0"; r[x == "H"] <- "0/1";  r[x == "B"] <- "1/1";  r[is.na(x)] <- "./.";  return(r)
}

fromVCFcode.geno <- function(x){
  r <- rep(NA,length(x));  r[x == "0|0"] <- "AA";  r[x == "0|1"] <- "AB"; r[x == "1|0"] <- "BA";  r[x == "1|1"] <- "BB";  r[x == ".|."] <- NA;  return(r)
}

fromVCFcode.AHB <- function(x){
  r <- rep(NA,length(x));  r[x == "0|0"] <- "A";  r[x == "0|1"] <- "H"; r[x == "1|0"] <- "H";  r[x == "1|1"] <- "B";  r[x == ".|."] <- NA;  return(r)
}

fromVCFcode.AHBp <- function(x){
  r <- rep(NA,length(x));  r[x == "0|0"] <- "A";  r[x == "0|1"] <- "H0"; r[x == "1|0"] <- "H1";  r[x == "1|1"] <- "B";  r[x == ".|."] <- NA;  return(r)
}

fromVCFcode.AHBn <- function(x){
  r <- rep(NA,length(x));  r[x == "0|0"] <- NA;  r[x == "0|1"] <- "H0"; r[x == "1|0"] <- "H1";  r[x == "1|1"] <- NA;  r[x == ".|."] <- NA;  return(r)
}
