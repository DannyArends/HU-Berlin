
setwd("D:/Edrive/Goat/DNA/SihamAnalysis")
tiff("Combined3Figs.tiff", width = 8000, height = 4000, compression  = "lzw", res=300)
  layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE),widths=c(3,1))
  par(cex=1.0)
  plot(dendrogram2.col, main = "Nei's genetic distance",cex.axis=1.2, cex.main=1.4, las=2)
  aa <- plotStructure(stmatrix, popinfo, TRUE, FALSE)

  plot(c(-50,100), c(-100,150), col = cols[as.character(viewn[as.character(groups)])],pch = 19, xlab=paste0("PC1 ",pca1), ylab=paste0("PC2 ",pca2), 
      t = 'n',xaxt='n', yaxt='n', main="Principal component analysis", cex.lab=1)
  axis(1, at = seq(-50, 100, 20),cex.axis=1.2)
  #abline(v = seq(-50, 100, 15), col="gray", lty=2)
  axis(2, at = seq(-100, 150, 20),las=2,cex.axis=1.2)
  #abline(h = seq(-100, 150, 50), col="gray", lty=2)
  points(pcares$x[,1], pcares$x[,2], col = cols[as.character(viewn[as.character(groups)])], pch = types[viewn[as.character(groups)]], cex=1.0)
  legend("topright", c("Taggar", "Desert", "Nilotic", "Nubian"), col=cols, pch=types, bg="white", cex=1.2)
dev.off()
