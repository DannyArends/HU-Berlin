# Analysis of Atlas Data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#
# Example for Stefan Kärst: Change 0 by NA in a matrix

stefan <- matrix(sample(0:3,1000,replace=TRUE),100,10)
rownames(stefan) <- paste("r0wnames",1:100)
colnames(stefan) <- paste("c0lnames",1:10)
apply(stefan,2,function(x){ x[which(x==0)] <- NA; return(x) })