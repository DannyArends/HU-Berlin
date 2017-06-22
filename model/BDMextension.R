mcM <- NULL
mcMo <- NULL

for(x in 1:1000){
  N <- 1000

  m1 <- sample(c("AA", "AT", "AT", "TT"),N, replace=TRUE)
  m2 <- sample(c("CC", "CG", "CG", "GG"),N, replace=TRUE)

  cM <- paste0(m1,":",m2)
  cMo <- cM
  op <- par(mfrow=c(1,2))
  plot(as.factor(cM))
  table(as.factor(cM)) / N

  i <- 1
  for(e in cM){
    if(e == "AA:CC" || e == "AA:CG"){
      new <- paste0(sample(c("AA", "AT", "AT", "TT"), 1, replace=TRUE), ":", sample(c("CC", "CG", "CG", "GG"), 1, replace=TRUE))
      while(new %in% c("AA:CC","AA:CG")){
        new <- paste0(sample(c("AA", "AT", "AT", "TT"), 1, replace=TRUE), ":", sample(c("CC", "CG", "CG", "GG"), 1, replace=TRUE))
      }
      #cat(i,"",e, "->",new,"\n")
      cM[i] <- new
    }
    i <- i + 1
  }

  plot(as.factor(cM))
  table(as.factor(cM)) / N

  tcM <- table(cM) / N
  tcMo <- table(cMo) / N

  mcM <- rbind(mcM, tcM)
  mcMo <- rbind(mcMo, tcMo)
  
  for(e in names(tcM)){
     cat(e, tcM[e] - tcMo[e], "\n")
  }
}
tcM <- apply(mcM, 2, mean)
tcM.sd <- apply(mcM, 2, sd)
tcMo <- apply(mcMo, 2, mean)
tcMo.sd <- apply(mcMo, 2, sd)
  
op <- par(mfrow=c(1,1))
op <- par(cex=1.4)

L <- 2* length(names(tcMo))
plot(c(0.5, L+0.5), c(0, 5/16 * 110), t = "n", xaxt='n', yaxt='n', xaxs="i", yaxs="i", xlab="Genotype", ylab="")
abline(v = seq(2.5, L, 2))
abline(h = seq(1/16 * 100, 40 , 1/16 * 100))
axis(1, at=seq(1.5,L,2), names(tcMo))
axis(2, at=seq(1/16 * 100, 5/16 * 100 , 1/16 * 100), c("1/16","1/8","3/16","1/4","5/16"), las=2)
i <- 1.0
for(e in names(tcMo)){
   rect(i,0,i+0.5,tcMo[e] * 100, col="orange")
   lines(c(i+0.25,i+0.25), c((tcMo[e] - tcMo.sd[e]) * 100, (tcMo[e] + tcMo.sd[e]) * 100))
   i <- i + 0.5
   rect(i,0,i+0.5,tcM[e] * 100, col="cornflowerblue")
   lines(c(i+0.25,i+0.25), c((tcM[e] - tcM.sd[e]) * 100, (tcM[e] + tcM.sd[e]) * 100))
   i <- i + 1.5
}
