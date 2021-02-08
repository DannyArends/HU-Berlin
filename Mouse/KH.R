# Redo figures of the KH

day1 <- rbind(c( 0, 2),
              c(20, 8),
              c(31,23),
              c( 6, 6),
              c(11,18),
              c(28,38),
              c( 2, 4),
              c( 0, 0))
day1.sem <- rbind(c( 0, 3),
                  c(14, 3),
                  c( 1, 7),
                  c( 2, 7),
                  c(11, 5),
                  c(28, 9),
                  c( 2, 4),
                  c( 0, 0))

plot(c(1,8), c(0,60), t = 'n', xaxt='n', las=2, xlab="Morphogenesis Stage", ylab="% of hair follicles in stage", yaxs="i")
abline(h = c(10, 20, 30, 40, 50), lty = 2)
points(day1[,1], col="blue", t = 'l',lwd=3)
for(x in 1:nrow(day1.sem)){

}
points(day1[,2], col="gold", t = 'l', lwd=3)

axis(1, at = 1:8, 1:8, tick = FALSE)
axis(1, at = 1:8+0.5, 1:8, labels = FALSE)