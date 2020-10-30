
toLOD <- function(models, term = "gt"){
  anovas <- lapply(models, anova)
  pvals <- unlist(lapply(anovas, function(x){ return(x[term, "Pr(>F)"]) }))
  lods <- -log10(pvals)
  return(lods)
}

chrLength <- function(map){
  lengths <- c()
  for(chr in chroms){ lengths <- c(lengths, max(map[map[,1] == chr,2], na.rm=TRUE)); }
  names(lengths) <- chroms
  return(lengths)
}

addCumLength <- function(map, chrLengths, gap = 50000000){
  map <- cbind(map, cumPos = NA)
  map[which(map[,1] == chroms[1]),"cumPos"] <- map[which(map[,1] == chroms[1]), 2]
  for(x in 2:length(chroms)){
    map[which(map[,1] == chroms[x]), "cumPos"] <- as.numeric(map[which(map[,1] == chroms[x]),2]) + as.numeric(sum(chrLengths[1:(x-1)])) + as.numeric(((x-1) * gap))
  }
  return(map)
}

chrCenter <- function(cmax, cmin){ return(((cmax - cmin) / 2) + cmin) }

chrCenters <- function(map, chrLengths, gap = 50000000){
  centers <- c()
  centers <- c(centers, chrCenter(max(map[which(map[,1] == chroms[1]), 2]), min(map[which(map[,1] == chroms[1]), 2])))
  for(x in 2:length(chroms)){
    v <- as.numeric(map[which(map[,1] == chroms[x]),2]) + as.numeric(sum(chrLengths[1:(x-1)])) + as.numeric(((x-1) * gap))
    centers <- c(centers, chrCenter(max(v), min(v)))
  }
  return(centers)
}