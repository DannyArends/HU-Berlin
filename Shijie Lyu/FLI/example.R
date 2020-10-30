pvals <- runif(10000, 0, 0.9)

rankedData <- rank(pvals, ties.method="first") + 0.5
plot(-log10(rankedData / (max(rankedData) + 1)), -log10(pvals))
abline(a=0, b=1)


lambdaB <- round(median(qchisq(pvals, 1, lower.tail = FALSE), na.rm=TRUE) /  qchisq(0.5, 1),3)
chiSq <- qchisq(pvals, 1, lower.tail = FALSE)
p.deflated1 <- pchisq(chiSq / lambdaB, 1, lower.tail = FALSE)
lambdaA1 <- round(median(qchisq(p.deflated1, 1, lower.tail = FALSE),na.rm=TRUE) /  qchisq(0.5, 1),3)

cat("Before:", lambdaB, "After", lambdaA1, "\n")
plot(pvals, p.deflated1, cex=2)

chi <- qnorm(1 - (pvals / 2)) 
p.deflated2 = pchisq((chi^2)/lambdaB, df = 1, lower = FALSE)
lambdaA2 <- round(median(qchisq(1.0 - p.deflated2, 1),na.rm=TRUE) /  qchisq(0.5, 1),3)

cat("Before:", lambdaB, "After", lambdaA2, "\n")
points(pvals, p.deflated1, col='red')