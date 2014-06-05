# simplePlot.R
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified May, 2014
# first written May, 2014
#
# Basic line plot, for cattle data

setwd("D:/safa")
cattlenumbers <- read.csv("cattlenumber.csv")

plot(cattlenumbers, type = 'l', col = "red", lwd = 3)
