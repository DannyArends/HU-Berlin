setwd("C:/Users/Arends/Downloads")

mdata <- read.csv("COVID-19_aantallen_gemeente_cumulatief(2).csv", sep = ";")

groningen <- mdata[which(mdata[,2] == "GM0014"),]
utrecht <- mdata[which(mdata[,2] == "GM0344"),]

dates <- groningen[-1,1]
dates <- gsub("2020-", "", unlist(lapply(strsplit(dates, " "), "[",1)))
dates <- dates[(length(dates) - 30): length(dates)]

plot(c(1,length(dates)), c(-1,25), t = 'n', xaxt='n', las=2)
axis(1, at = 1:length(dates), dates, las=2)
abline(v=15, col = "darkgreen")
abline(v=22, col = "red")

gr <- unique(mdata[which(mdata[, "Province"] == "Groningen"),2])
ut <- unique(mdata[which(mdata[, "Province"] == "Limburg"),2])

mm <- c()
for(x in gr) {
  if(x != ""){
    groningen <- mdata[which(mdata[,2] == x),]
    groningen.smooth <- smooth(diff(groningen[, "Total_reported"]))
    groningen.smooth <- groningen.smooth / mean(groningen.smooth)
    groningen.smooth <- groningen.smooth[(length(groningen.smooth) - 30): length(groningen.smooth)]
    mm <- rbind(mm, groningen.smooth)
    points(groningen.smooth, t = 'l', col = "green", lty=2)
  }
}
points(apply(mm,2, mean), t = 'l', lwd=3, col = "darkgreen")

mm <- c()
for(x in ut){
  if(x != ""){
    utrecht <- mdata[which(mdata[,2] == x),]
    utrecht.smooth <- smooth(diff(utrecht[, "Total_reported"]))
    utrecht.smooth <- utrecht.smooth / mean(utrecht.smooth)
    utrecht.smooth <- utrecht.smooth[(length(utrecht.smooth) - 30): length(utrecht.smooth)]
    mm <- rbind(mm, utrecht.smooth)
    points(utrecht.smooth, t = 'l', col="red", lty=2)
  }
}
points(apply(mm,2, mean), t = 'l', lwd=3, col = "red")
