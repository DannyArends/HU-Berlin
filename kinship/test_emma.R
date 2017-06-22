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
    return(list(flags,snps, K))
}

SNPs <- matrix(0,10,3)
SNPs[1:5,1] <- 0.5
emma.kinship2(SNPs, "additive")

SNPs[1:5,3] <- 1
emma.kinship2(SNPs, "additive")