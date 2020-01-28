
# Continue with the numeric genotypes on the server
setwd("~/AfricanGoats")
library(parallel)

data.num <- read.csv(file = "Africa.Swiss.Adaptmap.merged.num.matrix.txt", sep = "\t")

calcEuclidean <- function(x){
  res <- rep(NA, ncol(data.num))
  for(y in 1:ncol(data.num)){
    res[y] <- sqrt(sum((data.num[,x] - data.num[,y])^2, na.rm = TRUE))
  }
  return(res)
}
cl <- makeCluster(getOption("cl.cores", 25)) # Use 25 cores on the server
clusterExport(cl, "data.num")
res <- parLapply(cl, 1:ncol(data.num), calcEuclidean)
stopCluster(cl)

dist.mat <- matrix(NA, ncol(data.num), ncol(data.num), dimnames =list(colnames(data.num), colnames(data.num)))
for(r in 1:length(res)){
  dist.mat[r,] <- res[[r]]
}

write.table(dist.mat, file = "Africa.Swiss.Adaptmap.distancematrix.txt", sep = "\t", quote=FALSE)
