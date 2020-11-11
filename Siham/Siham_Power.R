
is.scalar <- function(v) { is.numeric(v) && length(v)==1; }
# Accessory function: Check whether the input is numeric scalar.
# Used by other functions for power calculation in checking input consistency.
# The user doesn't need to use this function directly.


power_n_hsq <- function(n, qsq, pval = 5E-8) {
# Calculates power based on a range of values of n and qsq. pval has to be fixed.
# Parameters as defined above.
# n and hsq can be vectors.
# pval a single numeric value. Default is the common genome-wide significance threshold 5E-8.
# Returns pow, a matrix containing power values, with n along rows and qsq along columns.
# Can contain NA if NCP is negative, NA or infinite, which may happen e.g. if qsq=1.
# Example call: pow = power_n_hsq(n = (1:5)*1000, qsq = (1:10)/100, pval=5E-8)

if (missing(n)) { stop("Parameter n not found."); }
if (missing(qsq)) { stop("Parameter qsq not found."); }
if ( !( is.vector(n) && is.numeric(n) ) ) { stop("Parameter n not a numeric vector."); }
if ( !( is.vector(qsq) && is.numeric(qsq) ) ) { stop("Parameter qsq not a numeric vector."); }
if ( !( is.scalar(pval) ) ) { stop("Parameter pval not a numeric scalar."); }
if ( any( n<0 | !is.finite(n) ) ) { stop("Parameter n has unacceptable values."); }
if ( any( qsq<0 | qsq>1 | !is.finite(qsq) ) ) { stop("Parameter qsq has unacceptable values."); }
if ( pval<=0 | pval>=1 | !is.finite(pval) ) { stop("Parameter pval has unacceptable value."); }

th=qchisq(pval,df=1,lower.tail=F); # Significance threshold for chi-square, corresponding to P-value threshold

pow=matrix(NA,nrow=length(n),ncol=length(qsq)); # Pre-allocalte matrix to store power values
rownames(pow)=n; colnames(pow)=qsq; # Assign row and column names

for (i in 1:length(n)) {
	for (j in 1:length(qsq)) {
		ncp=n[i]*qsq[j]/(1-qsq[j]); # Calculate NCP parameter
		if (is.finite(ncp) && ncp>=0) { pow[i,j]=pchisq(th,df=1,lower.tail=F,ncp=ncp); } # Calculate power when NCP is finite and >=0
}}

return(pow); # Return calculated power matrix
}

aa <- power_n_hsq(seq(25, 1000, 25), (1:10)/100 ,2e-06)
colfunc <- colorRampPalette(c("red","gold","springgreen","royalblue"))
colz <- colfunc(ncol(aa))

plot(c(0,1000), y = c(0, 1), xlab="Sample Size", ylab = "Power", main= "Power analysis GWAS (threshold = 2e-06)", sub = "Computation based on DOI: 10.1016/j.ajhg.2017.06.005", t = 'n', las=2,xaxt='n')
axis(1, at = seq(0,1000,50), seq(0,1000,50))
i <- 1
apply(aa, 2, function(x){
  points(rownames(aa), x, t = 'b', col = colz[i], pch=19)
  i <<- i + 1
})
abline(h = 0.8, lty=2, col = "gray")
abline(v = 500, lty=3, col = "gray")
legend("topleft", paste0(round(as.numeric(colnames(aa)) * 100), "%"), lwd=1, pch=19, col=colz, bg="white", title = "Variance explained by SNP", ncol = 2)

