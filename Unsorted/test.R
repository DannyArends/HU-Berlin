
pvals <- runif(10000) * 0.002

lambdaB <- round(median(qchisq(pvals, 1, lower.tail = FALSE), na.rm=TRUE) /  qchisq(0.5, 1),3)
chiSq <- qchisq(pvals, 1, lower.tail = FALSE)
p.deflated1 <- pchisq(chiSq / lambdaB, 1, lower.tail = FALSE)
lambdaA1 <- round(median(qchisq(p.deflated1, 1, lower.tail = FALSE),na.rm=TRUE) /  qchisq(0.5, 1),3)

plot(pvals, p.deflated1 )