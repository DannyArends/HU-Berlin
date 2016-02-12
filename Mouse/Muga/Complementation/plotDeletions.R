setwd("E:/Mouse/ClassicalPhenotypes/Complementation/RawData")

df <- read.table("dfraw.txt", header =T, sep = "\t", colClasses="character")
bbs7 <- read.table(file="Bbs7_visual.txt", sep="\t", na.strings = "NA", stringsAsFactors=F, header=TRUE ) 

plot.start <- 36598000
plot.end <- 36604000

d <- data.frame(seq(plot.start, plot.end, by=1000))
colnames(d) <- "POS"
d$level <- 0

#dev.new(width=1600, height=800)

plot(d$POS, d$level, ylim = c(-1.5,2.5), pch=".", col= "white", xaxt ="n", main="Bbs7 - Exon/Intron structure", xlab="Position on Chr3", yaxt = "n", ylab="")#, yaxt ="n", ylab = "") 
lines <- seq(plot.start,plot.end,100)
abline(v=lines, col= "grey90")
lines <- seq(plot.start,plot.end,1000)
abline(v=lines, col= "grey75")
#mtext("only SNPs that diff. between Leans & BFMI", side = 3, cex = 0.7)
#mtext("BFMI GT from raw data (low Qual)", side = 3, cex = 0.7)
axis(1, d$POS, labels = prettyNum(d$POS, big.mark = ","))
segments(36571000, 1, 36614500, 1, col = "steelblue", lty = "dotted", lwd = 2)    # rand
segments(36573142, 1, 36613477, 1, col = "steelblue", lty = 1, lwd = 2)   # Bbs7 max
rect(36572000, 6.4, 36614700, 10,      density= NA, col = "grey60", border = NULL, angle = 135, lwd = 0.1)
rect(36572000, 3   , 36614700, 6.3 , density= NA, col = "orange", border = NULL, angle = 45, lwd = 0.1)

## show intron/exons
is.odd <- function(x) x %% 2 != 0 
for(x in 1:nrow(bbs7)){
  if(is.odd(x)){
    segments(bbs7[x,"Start"],1, bbs7[x,"End"],1, col = "orange2", lty = 1, lwd = 10, lend="butt") # (x/10+0.1) Introns
    text(mean(c(bbs7[x,"Start"],bbs7[x,"End"])),(x%%3)*0.4-0.1,bbs7[x,"No"], cex=0.8, font=2)
    text(mean(c(bbs7[x,"Start"],bbs7[x,"End"])),(x%%3)*0.4-0.3,bbs7[x,"Length"], cex=0.8)  
    segments(mean(c(bbs7[x,"Start"],bbs7[x,"End"])),0.9, mean(c(bbs7[x,"Start"],bbs7[x,"End"])), (x%%3) * 0.4,  col = "grey30", lty = "dotted", lwd = 2, lend="butt")
  } else {
    segments(bbs7[x,"Start"],1, bbs7[x,"End"],1, col = "red", lty = 1, lwd = 20, lend="butt") # Exons
    text(mean(c(bbs7[x,"Start"],bbs7[x,"End"])),(x%%3)*0.4+1.6,bbs7[x,"No"], cex=0.8, font=2)
    text(mean(c(bbs7[x,"Start"],bbs7[x,"End"])),(x%%3)*0.4+1.45,bbs7[x,"Length"], cex=0.8)
    segments(mean(c(bbs7[x,"Start"],bbs7[x,"End"])),1.2, mean(c(bbs7[x,"Start"],bbs7[x,"End"])), (x%%3)*0.4+1.3,  col = "grey30", lty = "dotted", lwd = 2, lend="butt")
  }
}

#text(mean(d$POS), -1, "Bbs7", font = 2)
segments(36599687, 1, 36601265, 1, col = "grey80", lty = par("lty"), lwd = 10, lend="butt")       # Bbs7 del Intron8
segments(36599798, 1, 36600147, 1, col = "steelblue", lty = par("lty"), lwd = 10, lend="butt")    # Bbs7 del CTCF
text(mean(c(36599798,36600147)), 1.3, "CTCF", col= "steelblue", cex=0.8)
segments(36613669, 1, 36613478, 1, col = "gold", lty = par("lty"), lwd = 10, lend="butt")         # 5'
segments(36573143, 1, 36572792, 1, col = "gold", lty = par("lty"), lwd = 10, lend="butt")         # 3'
segments(36613563, 1, 36613583, 1, col = "grey30", lty = par("lty"), lwd = 10, lend="butt")       # Bbs7 del prom, too small to show , in echt: 36613573!!!
text(mean(c(36613563,36613573))+100, 1.3, "del", col= "gray30", cex=0.8) 

# Lines
segments(36600373, -0.3, 36601300, -0.3, col = "orange", lty = "dotted", lwd = 2, lend="butt")
segments(36600026, -0.6, 36601300, -0.6, col = "orange", lty = "dotted", lwd = 2, lend="butt")
segments(36599493, -0.9, 36601300, -0.9, col = "orange", lty = "dotted", lwd = 2, lend="butt")
segments(36599493, -1.2, 36601300, -1.2, col = "orange", lty = "dotted", lwd = 2, lend="butt")

# Sequenzierung
segments(36601267, -0.3, 36601389, -0.3, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.3, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1001_vorDelup2:174-296
segments(36600373, -0.3, 36600515, -0.3, col = "orange", lty = par("lty"), lwd = 4, lend="butt")#;text(36601389+10,-0.3, "M1001", srt=0, adj = 0, cex=0.8) # sequenzierung  M1001_vorDelup:300-441
segments(36601267, -0.3, 36601389, -0.3, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.3, "NZO.M1001", srt=0, adj = 0, cex=0.8)  # sequenzierung M1001_vorDelup2:174-296

segments(36601267, -0.6, 36601389, -0.6, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1000_vorDelup2:173-295
segments(36600370, -0.6, 36600515, -0.6, col = "orange", lty = par("lty"), lwd = 4, lend="butt")#;text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1000_vorDelup:300-444
segments(36601267, -0.6, 36601389, -0.6, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8)      # sequenzierung M1000_vorDelup2:173-295
segments(36599525, -0.6, 36600114, -0.6, col = "darkorange3", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8) # sequenzierung NZO1000_hiDELlow:5-597
segments(36600026, -0.6, 36600178, -0.6, col = "darkorange3", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8) # sequenzierung NZO1000_DELlow:4-156
segments(36600242, -0.6, 36600516, -0.6, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8)      # sequenzierung NZO1000_prDEL2up:30-303
segments(36600242, -0.6, 36600515, -0.6, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8)      # sequenzierung NZO1000_vorDELup:295-567
segments(36601390, -0.6, 36601519, -0.6, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601519+10,-0.6, "NZO.M1000", srt=0, adj = 0, cex=0.8) # sequenzierung NZO1000_vorDELup:21-149
segments(36601267, -0.6, 36601389, -0.6, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601389+10,-0.6, "", srt=0, adj = 0, cex=0.8)      # sequenzierung NZO1000_vorDELup:172-294

segments(36601265, -0.9, 36601544, -0.9, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601544+10,-0.9, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1544_vorDelup:1-277
segments(36599496, -0.9, 36599686, -0.9, col = "orange", lty = par("lty"), lwd = 4, lend="butt")#;text(36601544+10,-0.9, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1544_vorDelup:280-470
segments(36599601, -0.9, 36599686, -0.9, col = "orange", lty = par("lty"), lwd = 4, lend="butt")#;text(36601544+10,-0.9, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1544_vorDelup2:280-365
segments(36601265, -0.9, 36601544, -0.9, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601544+10,-0.9, "BFMI.M1544", srt=0, adj = 0, cex=0.8)  # sequenzierung M1544_vorDelup:1-277
segments(36599493, -0.9, 36599686, -0.9, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36601544+10,-0.9, "", srt=0, adj = 0, cex=0.8)      # sequenzierung M1544_prDEL2up:44-237  
segments(36599526, -0.9, 36599686, -0.9, col = "darkorange3", lty = par("lty"), lwd = 4, lend="butt");text(36601544+10,-0.9, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1544_hiDELlow:6-165
segments(36601265, -0.9, 36601331, -0.9, col = "darkorange3", lty = par("lty"), lwd = 4, lend="butt");text(36601544+10,-0.9, "", srt=0, adj = 0, cex=0.8) # sequenzierung M1544_hiDELlow:168-234

segments(36599493, -1.2, 36599686, -1.2, col = "orange", lty = par("lty"), lwd = 4, lend="butt");text(36599686+10,-1.2, "", srt=0, adj = 0, cex=0.8)      # sequenzierung AKR1009_prDEL2up:45-238
segments(36599524, -1.2, 36599686, -1.2, col = "darkorange3", lty = par("lty"), lwd = 4, lend="butt");text(36599686+10,-1.2, "", srt=0, adj = 0, cex=0.8) # sequenzierung AKR1009_hiDELlow:1-162
segments(36601265, -1.2, 36601331, -1.2, col = "darkorange3", lty = par("lty"), lwd = 4, lend="butt");text(36601331+10,-1.2, "AKR1009", srt=0, adj = 0, cex=0.8) # sequenzierung AKR1009_hiDELlow:165-231
