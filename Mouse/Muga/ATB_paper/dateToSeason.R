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

phenonames <- c("mri42d_fatDlean", "mri56d_fatDlean", "mri70d_fatDlean",                                        # Fat / Lean 
                "d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71",                                  # Weights
                "GF1", "GF2", "total.GF",                                                                       # Gonadal fat
                "RF1", "RF2", "total.RF", "IF", "Muskel", "Leber", "BAT", "LD",                                 # Other phenotypes
                "mri42d_fat", "mri42d_lean", "mri42d_3", "mri42d_4",                                            # MRI day 42
                "mri56d_fat", "mri56d_lean", "mri56d_3", "mri56d_4",                                            # MRI day 56
                "mri70d_fat", "mri70d_lean", "mri70d_3", "mri70d_4")                                            # MRI day 70