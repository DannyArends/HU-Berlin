# \file Answers - Computing Fundamentals.R
#
# Copyright (c) 2014-2016, Danny Arends
# Last modified:		May, 2014
# Created:		      May, 2014
# 
# A copy of the GNU General Public License, version 3, is available
# at http://www.r-project.org/Licenses/GPL-3
# 
# Contains: Answers to the Computing Fundamentals

setwd("~/practicals/Computing Fundamentals") # Set working dir

# 0)
vHello <- "Hello World"

cat(vHello)         # Note, no carriage return
print(vHello)
vHello

# 1a)
unknown <- 4
if(unknown >= 0) print("Higher then or equal to 0")

# 1b)
unknown <- 4
if(unknown < 0 || unknown > 10) stop("Variable unknown is not between 0 and 10")

# 2)
forsum <- 0
for(x in 1:1000){
  forsum <- forsum + x
}

whilesum <- 0
x <- 1
while(x <= 1000){
  whilesum <- whilesum+x
  x <- x + 1
}

if(forsum == whilesum) print("They match, the outcome is: ", forsum)

# 3)

randomnr <- runif(1, 10)
cat(randomnr, "is ")
if(randomnr < 5) cat("lower")
if(randomnr > 5) cat("higher")
cat(" then 5")

# 4)

for(x in 1:12){
  n <- 1
  while(n < x){
    cat("#")
    n = n + 1
  }
  cat("\n")
}

# 5)
uniformMatrix <- matrix(runif(100), nrow=20, ncol=5)
gaussianMatrix <- matrix(rnorm(100), nrow=20, ncol=5)

histogram(uniformMatrix[,1])                            # Generate plot
histogram(gaussianMatrix[,1])                           # Generate plot

write.table(uniformMatrix, file="uniformMatrix.txt")    # Write to file
write.table(gaussianMatrix, file="gaussianMatrix.txt")  # Write to file

# 6)
cat("I say: \"Escaping stuff is \'great\', but \\ and / might be a nuisance.\n", file="escape.txt")
cat("You are correct, sir!", file="escape.txt", append=TRUE)


# Additional

# AA1) This is actually a question people use in elementary school to test children’s ability for mathematics... :)
# Use the following, imagine that we take the row two times, reverse one of the rows and put it below the other one.
# Then add them up though the columns
# 1 + 2 + 3 + 4
# 4 + 3 + 2 + 1
# ============== +
# 5 + 5 + 5 + 5
# We end up with 4 times 5, or x * (x+1), however, we took the row two times so we need to divide the answer by two
# so 5 * 4 / 2, or (x * (x+1)) / 2

x <- 1000
smartsum = (x * (x+1)) / 2

# AA2)
as.integer(as.factor(c(3, 4, 5)))

# AA3) For a 24 x 24 grid

w <- 0                                        # Width of the current row
for(x in 1:24){                               # Current x-coordinate
  y <- 1                                      # Current y-coordinate
  while(y < 24){
    if(y >= (12 - w) && y <= (12 + w)){       # If we are between the middle - w and middle + w
      cat("#")                                # Print a #
    }else{
      cat(" ")                                # Print a space
    }
    y = y + 1
  }
  if(x < 12) w = w + 1                        # If we're below 12, increase our width
  if(x > 12) w = w - 1                        # If we're above 12, decrease our width
  cat("\n")                                   # End the current line
}

# To make it more general, replace 24 with size, and 12 by size/2, now we can print as small / big as we want
