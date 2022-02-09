### TRD figure 2, version from Paula (18.03.2021)
# load libraries and set wd
options(width=180, stringsAsFactors=F)
library(ggplot2)
library(reshape2)
setwd("D:/Edrive/Mouse/DNA/MegaMuga/TRDredo")

# load chiscores, transform to LODscores and keep only signficant ones
chiscores    = read.table("ChiSqScores_Incompatible.txt", header=T, check.names=F)
LODscores    = apply(chiscores, 1, function(x){-log10(pchisq(x, 1, lower.tail=FALSE))})
threshold = -log10(0.05 / (1008  * 1008  * 0.5))
LODscores[which(LODscores < threshold, arr.ind=T)]= NA
LODscores = melt(LODscores, na.rm=T)  # "melting"  data.frame in order to be suitable for ggplot

# load alleles and add them to the LODscores data.frame
d = read.table("AlleleCombo_Incompatible.txt", header=T, check.names=F) # 62 62
LODscores$A1 = NA
LODscores$A2 = NA
for (i in 1:nrow(LODscores)){
	alleles = d[LODscores$Var1[i], LODscores$Var2[i]]
	LODscores$A1[i] = gsub("\\|.+", "", alleles)
	LODscores$A2[i] = gsub(".+\\|", "", alleles)
}


# re-factor genomic positions (as only 29 out of 62 are used) and make them numeric to get x- any y-coordinates for the plot 
LODscores$Var1 = factor(LODscores$Var1, intersect(levels(LODscores$Var1), unique(LODscores$Var1)))
LODscores$Var2 = factor(LODscores$Var2, intersect(levels(LODscores$Var2), unique(LODscores$Var2)))
LODscores$Var1num = as.numeric(LODscores$Var1)
LODscores$Var2num = as.numeric(LODscores$Var2)
# melting data.frame again, as a rectangle has to be made for each allele separately
LODscores2 = melt(LODscores, id.vars= setdiff(colnames(LODscores), c("A1", "A2")), variable.name="allele", value.name="allelename")
# adjusting x- and y-cooorinates to be shifted 0.25 left if first allele or 0.25 right if second allele
LODscores2$x = LODscores2$Var1num
LODscores2$x[LODscores2$allele=="A1"] = LODscores2$x[LODscores2$allele=="A1"] - 0.25
LODscores2$x[LODscores2$allele=="A2"] = LODscores2$x[LODscores2$allele=="A2"] + 0.25
# factor allele to get right order in plot legend
LODscores2$allelename = factor(LODscores2$allelename, levels=c("B6N", "H", "BFMI"))

# generate data.frame for empty squares (estimating mean x-/y-cooridnate per chr and size of square )
emptysquares = data.frame(name= levels(LODscores$Var1), order= 1:29, chr = gsub(":.+", "", levels(LODscores$Var1)))
number       = as.data.frame(table(emptysquares$chr), stringsAsFactors=F)      # size of square == frequency of chr
emptysquares = merge(emptysquares, number, by.x="chr", by.y="Var1")
emptysquares = aggregate(order ~ Freq + chr, emptysquares, mean)  # mean x-/y-cooridnate

## make plot (order of the single elements which you add to the plot plays a role)
ggplot(LODscores2, aes(Var1, Var2)) + theme_bw() +  # Load data 
	# add signficant hits
	geom_tile(aes(x= x, y=Var2num, fill=allelename), color="white") +
	# add vertical/horizontal grey lines in background
	geom_vline(xintercept= seq(0.5,28.5,1), color="white", size=0.1, linetype = 1)+
	geom_vline(xintercept= seq(0.5,28.5,1), color="grey50", size=0.1, linetype = 3)+
	geom_hline(yintercept= seq(0.5,28.5,1), color="white", size=0.1, linetype = 1)+
	geom_hline(yintercept= seq(0.5,28.5,1), color="grey50", size=0.1, linetype = 3)+
	# add emptysquares in grey
	geom_tile(data=emptysquares, aes(x=order, y=order, width=Freq, height=Freq), fill="grey")+
	# add vertical/horizontal lines between chromosomes
	geom_vline(xintercept= c(which(!duplicated(gsub(":.+","",levels(LODscores2$Var1)))),30)-0.5, size=0.5)+
	geom_hline(yintercept= c(which(!duplicated(gsub(":.+","",levels(LODscores2$Var1)))),30)-0.5, size=0.5)+
	# rotate x-axis labels by 90 degrees and remove default grid lines
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), panel.grid.minor = element_blank(),panel.grid.major = element_blank() )+
	# make thicks for x-/y-axis, add labels (genomic positions), set x-/y-axis limits, and remove space between the axis and the actual plot
	scale_x_continuous(breaks=1:29, labels= levels(LODscores$Var1), limits=c(0.5,29.5), expand=c(0,0))+
	scale_y_continuous(breaks=1:29, labels= levels(LODscores$Var1), limits=c(0.5,29.5), expand=c(0,0))+
	# set legend name, color for rectangles (order: B6N, H, BFMI), and fix coordinates so that the plot is quadratic
	scale_fill_manual(name="", values= toupper(c("cornflowerblue", "yellowgreen", "orange")))+coord_fixed(3/4)+
	# add x- and y- labels
	xlab("") + ylab("")

# save plot
ggsave("figure2.png", width=8, height=8, dpi =600)
