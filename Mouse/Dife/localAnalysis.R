setwd("E:/Mouse/DNA/Sequencing/DifeMouse/")

samples <- c("ext_L7254", "ext_L7256", "ext_L7257", "ext_L7255", "ext_L7258")
names(samples) <- c("NZO", "BFMI-S1", "BFMI-S2", "SJL", "BFMI-S12")

snpdata <- read.csv("out.txt", sep = "\t", header=TRUE)

