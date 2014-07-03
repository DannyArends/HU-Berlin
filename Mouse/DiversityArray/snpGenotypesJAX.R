# Analysis of JAX data from the mouse diversity array
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified June, 2014
# first written June, 2014
#

setwd("E:/Mouse/DNA/DiversityArray/")

calldataJAX    <- read.table("Analysis/MeasurementsJAX.txt", sep="\t", header=TRUE, colClasses=c("character"))                        # Load the unannotated data from JAX
calldataAtlas  <- read.table("Analysis/measurementsAtlas_annotated.txt", sep="\t", header=TRUE, colClasses=c("character"))            # Load the annotated Atlas dataset

JaxInAtlas     <- which(calldataJAX[,"JAX_ID"] %in% calldataAtlas[,"JAX_ID"])                                                         # Which JAX measurements are in Atlas
cat("We can merge:", length(JaxInAtlas),"/", nrow(calldataJAX),"probes from", nrow(calldataAtlas),"probes from JAX to Atlas\n")

calldataJAX     <- calldataJAX[JaxInAtlas, ]                                                                                          # Throw away 80.000 measurements from JAX

JAXinAtlas  <- match(calldataAtlas[,"JAX_ID"], calldataJAX[,"JAX_ID"])                                                                # Align JAX with the Atlas data
calldataJAX <- calldataJAX[JAXinAtlas,]
write.table(cbind(calldataAtlas[,1:9], calldataJAX[,-1]), file="Analysis/measurementsJAX_annotated.txt", sep="\t", row.names=FALSE)   # Write out the annotated numeric JAX genotypes

calldata <- cbind(calldataAtlas, calldataJAX[,-1])                                                                                    # Combine JAX and Atlas
write.table(calldata, file="Analysis/measurementsALL_annotated.txt", sep="\t", row.names=FALSE)                                       # Write out the annotated numeric combined genotypes