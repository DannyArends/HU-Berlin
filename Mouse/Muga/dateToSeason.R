# Convert dates to seasons
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Aug, 2014
# first written Aug, 2014

getSeason <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES), ".", fixed=TRUE), "[", 2)))
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5]                  <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8]                  <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11]                 <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2]  <- "Winter"
  return(ret)
}
