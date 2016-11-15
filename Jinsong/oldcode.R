

# OLDER CODE
get.eff.lme(phe.sex.age.adj, phe.name, "Strain")

lvls <- seq(0,1,length.out = nlevels(as.factor(phenotypes[,"Sex"])))
cols <- rgb(lvls, rep(1,length(lvls)), 1.0 - lvls, 0.5)
phe.adj.sex <- get.residuals(get.eff(phenotypes, phe.name, "Sex")) +  mean(as.numeric(phenotypes[,phe.name]), na.rm = TRUE)

op <- par(mfrow=c(2, 1))

boxplot(as.numeric(phenotypes[,phe.name]) ~ as.factor(phenotypes[,"Sex"]), col = cols, notch = TRUE)
boxplot(phe.adj.sex ~ as.factor(phenotypes[,"Sex"]), add=TRUE, col = cols, notch = TRUE)

phe.adj.logAge <- get.residuals(get.eff.num(phenotypes, phe.name, "LogAGE")) +  mean(as.numeric(phenotypes[,phe.name]), na.rm = TRUE)

plot(as.numeric(phenotypes[,phe.name]) ~ as.numeric(phenotypes[,"LogAGE"]), col = rgb(1.0, 0.0, 0.0, 0.7), pch = 18, ylab = phe.name, main = phe.name)
points(phe.adj.logAge ~ as.numeric(phenotypes[,"LogAGE"]), col = rgb(0.5, 1.0, 0.2, 0.7), pch = 18)
legend("topright", c("Before Age correction", "After Age correction") ,col =c(rgb(1.0, 0.0, 0.0, 0.7), rgb(0.5, 1.0, 0.2, 0.7)), pch = 18, bg  = "white")



