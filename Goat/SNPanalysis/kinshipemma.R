emma.kinship2 <- function(snps, method="additive", use="all") {
    n0 <- sum(snps==0,na.rm=TRUE)
    nh <- sum(snps==0.5,na.rm=TRUE)
    n1 <- sum(snps==1,na.rm=TRUE)
    nNA <- sum(is.na(snps))

    stopifnot(n0+nh+n1+nNA == nrow(snps) * ncol(snps))

    if ( method == "dominant" ) {
        flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
        snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
    }
    else if ( method == "recessive" ) {
        flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
        snps[!is.na(snps) && (snps == 0.5)] <- flags[!is.na(snps) && (snps == 0.5)]
    }
    else if ( ( method == "additive" ) && ( nh > 0 ) ) {
        dsnps <- snps
        rsnps <- snps
        flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) > 0.5),nrow(snps),ncol(snps))
        dsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]  # modified
        flags <- matrix(as.double(rowMeans(snps,na.rm=TRUE) < 0.5),nrow(snps),ncol(snps))
        rsnps[!is.na(snps) & (snps==0.5)] <- flags[!is.na(snps) & (snps==0.5)]  # modified
        snps <- rbind(dsnps,rsnps)
    }

    if ( use == "all" ) {
        mafs <- matrix(rowMeans(snps,na.rm=TRUE),nrow(snps),ncol(snps))
        snps[is.na(snps)] <- mafs[is.na(snps)]
    }
    else if ( use == "complete.obs" ) {
        snps <- snps[rowSums(is.na(snps))==0,]
    }

    n <- ncol(snps)
    K <- matrix(nrow=n,ncol=n)
    diag(K) <- 1

    for(i in 1:(n-1)) {
        for(j in (i+1):n) {
            x <- snps[,i]*snps[,j] + (1-snps[,i])*(1-snps[,j])
            K[i,j] <- sum(x,na.rm=TRUE)/sum(!is.na(x))
            K[j,i] <- K[i,j]
        }
    }
    return(K)
}

setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
numsnpdata <- read.csv("filtered_snps_numeric_NO_DN2.txt", sep="\t", check.names=FALSE)
samples    <- read.table("sampleinfo.txt", sep="\t")
locations  <- read.table("Sample_SNP_location_fixed.txt", sep="\t", header=TRUE, row.names=1)  # Phenotype data`
samples    <- cbind(samples, locations[rownames(samples),])
samples    <- cbind(samples, locationShort = as.character(unlist(lapply(strsplit(as.character(samples[,"Location"]), "_"), "[",1))))

asMAF <- ((numsnpdata - 1) / 2)
kinship <- emma.kinship2(asMAF)
colnames(kinship) <- colnames(numsnpdata)
rownames(kinship) <- colnames(numsnpdata)

breeds <- samples[,"Breed"]
names(breeds) <- rownames(samples)

colnames(kinship) <- breeds[colnames(kinship)]
rownames(kinship) <- breeds[rownames(kinship)]

colfunc <- colorRampPalette(c("green", "red", "white"))


heatmap(kinship, scale="none", col=colfunc(10))

snpdata <- read.table("filtered_snps.txt", sep="\t", check.names=FALSE, colClasses="character")
snpdata <- snpdata[,-which(colnames(snpdata) == "DN 2")] # Throw away the duplicate individual because it confuses STRUCTURE

res <- apply(snpdata,1,function(x){sum(is.na(x))})

snpdata <- snpdata[-which(res > 0), ]

slashc <- matrix(NA, nrow(snpdata), ncol(snpdata), dimnames=list(rownames(snpdata), colnames(snpdata)))
for(x in 1:nrow(snpdata)){
  slashc[x, ] <- unlist(lapply(strsplit(as.character(snpdata[x,]), ""), function(x){paste0(x[1], "/", x[2])}))
}
snpINFO <- snpinfo[rownames(snpdata),]

LDs <- vector("list", nrow(slashc))
names(LDs) <- rownames(slashc)
for(x in 1:nrow(slashc)){
  mChr <- as.character(snpINFO[x, "Chr"])
  mPos <- as.numeric(snpINFO[x, "Position"])
  nearby <- which(as.character(snpINFO[,"Chr"]) == mChr & as.numeric(snpINFO[, "Position"]) > (mPos - 1000000) &  as.numeric(snpINFO[, "Position"]) < (mPos + 1000000))
  mGeno <- genotype(slashc[x,])
  for(y in nearby){
     LDs[[x]] <- rbind(LDs[[x]], c(y, round(LD(mGeno,genotype(slashc[y,]))$"R^2",2)))
  }
  if(x %% 100 == 0) cat("Done",x,"/",nrow(slashc), "\n")
}


library(pegas)

library(poppr)
library(ape)

subNumMat <- apply(asMAF,2,as.numeric)


MAT <- t(subNumMat)
hasNA <- which(apply(MAT, 2, function(x){any(is.na(x))}))
MAT <- MAT[,-hasNA]
nloc <- ncol(MAT)


# Reynolds distance
denomi <- MAT %*% t(MAT)
vec <- apply(MAT, 1, function(x) sum(x * x))
D <- -2 * denomi + vec[col(denomi)] + vec[row(denomi)]
diag(D) <- 0
denomi <- 2 * nloc - 2 * denomi
diag(denomi) <- 1
D <- D/denomi
D <- sqrt(D)
colnames(D) <- breeds[colnames(D)]
rownames(D) <- breeds[rownames(D)]

D <- as.dist(D)

clusters <- hclust(D)
colnames(subNumMat) <- breeds[colnames(subNumMat)]
rownames(subNumMat) <- breeds[rownames(subNumMat)]

clusters2 <- hclust(dist(t(subNumMat),method ="manhattan"))

# Create colors
cols <- c("red", "blue", "orange", "black")
names(cols) <- c("T", "D", "Ni", "Nu")

viewn <- c("T", "D", "Ni", "Nu")
names(viewn) <- c("Tagg", "Dese", "Ni", "Nu")

clusters$labels <- viewn[clusters$labels]
clusters2$labels <- viewn[clusters2$labels]

labelCol <- function(x) {
  if (is.leaf(x)) {
    hclass <- as.character(attr(x, "label"))             # Fetch the class label
    hcol <- cols[hclass]                            # Determine color of the label
    cat(attr(x, "label"), hclass, hcol, "\n")
    attr(x, "nodePar") <- list(lab.col=hcol, cex=0.1)
  }
  return(x)
}

dendrogram1 <- as.dendrogram(clusters)
dendrogram2 <- as.dendrogram(clusters2)
dendrogram1.col <- dendrapply(dendrogram1, labelCol)
dendrogram2.col <- dendrapply(dendrogram2, labelCol)

                                                                  # Phylogram of all ALL lines
#plot(phylogram, type="fan")
#op <- par(mfrow=c(2,1))
plot(dendrogram1.col, main="Reynold distance")
plot(dendrogram2.col, main="Manhattan distance")
