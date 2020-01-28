
chromosomeOrder <- function(map){
  map <- cbind(map, pos = (map[, "Stop"] + map[, "Start"]) / 2)
  map <- map[sort(map[,"pos"], index.return = TRUE)$ix,]
  orderedmap <- c()
  for(chr in min(unique(map[, "chrN"])):max(unique(map[, "chrN"]))){
    orderedmap <- rbind(orderedmap, map[which(map[,"chrN"] == chr),c("chrN", "pos")])
  }
  colnames(orderedmap) <- c("chr", "pos")
  return(orderedmap)
}

tooCloselyRelated <- function(dist.mat, h = 10){
  diag(dist.mat) <- NA
  minDistances <- apply(dist.mat,1,min,na.rm=TRUE)
  indices <- which(minDistances < h)
  tooClose <- c()
  for(i in indices){
    tooClose <- rbind(tooClose, c(i, which.min(dist.mat[i,])))
  }
  torem <- c()
  for(x in 1:nrow(tooClose)){
    if(!(tooClose[x,1] %in% torem)){ torem <- c(torem, tooClose[x, 2]); }
  }
  return(torem)
}

# Function to color the dendrogram nodes based on the strain
labelCol <- function(x, LblAsBreed = FALSE){
  if(is.leaf(x)){
    label <- attr(x, "label") # Fetch label
    breed <- sample.annot[label, "Breed"]
    continent <- sample.annot[label, "Continent"]
    if(LblAsBreed) attr(x, "label") <- breed  # Set the label color based on the strain
    attr(x, "nodePar") <- list(lab.col = continents[continent])  # Set the label color based on the strain
  }
  return(x)
}

