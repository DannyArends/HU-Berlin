# Analysis of hetrozygousness in the 600k SNP chips
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Feb, 2015
# first written Feb, 2015

setwd("E:/Chicken/DNA/600KSNPChip/")

calldata <- read.table("GenotypeData.txt", na.strings="-1")
calldata <- calldata[which(rownames(calldata) %in% arrayAnnotation[,"Probe.Set.ID"]),]

arrayAnnotation <- read.table("Annotation/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", header=TRUE, colClasses="character", na.string="---")
arrayAnnotation <- arrayAnnotation[-which(is.na(arrayAnnotation[,"Physical.Position"])),c("Probe.Set.ID","Chromosome", "Physical.Position")]
arrayAnnotation <- arrayAnnotation[which(arrayAnnotation[,"Chromosome"] == 17), ]

testHetrozygousity <- function(calldata, arrayAnnotation, regionwidth = 100000){
  output <- NULL
  chromosomes <- as.character(unique(arrayAnnotation[,"Chromosome"]))
  for(chromosome in chromosomes){
    chrAnnotation <- arrayAnnotation[arrayAnnotation[,"Chromosome"] == chromosome,]
    maxlength <- max(chrAnnotation[,"Physical.Position"])
    regions <- cbind(start = seq(1,as.numeric(maxlength),regionwidth), end = seq(regionwidth,as.numeric(maxlength)+regionwidth,regionwidth))
    for(x in 1:nrow(regions)){
      snpInRegion <- chrAnnotation[which(as.numeric(chrAnnotation[,"Physical.Position"]) > regions[x,"start"] & as.numeric(chrAnnotation[,"Physical.Position"]) < regions[x,"end"]),"Probe.Set.ID"]
      if(length(snpInRegion) > 0){
        isHe <- sum(calldata[snpInRegion,] == 1,na.rm=TRUE)
        isHo <- sum(calldata[snpInRegion,] != 1,na.rm=TRUE)
      }else{
        isHe <- 0 ; isHo <- 1
      }
      output <- rbind(output, c(chromosome, regions[x,"start"], regions[x,"end"], round((isHe/isHo) * 100, d=2)))
    }
    cat("Done chromosome",chromosome,"\n")
  }
  colnames(output) <- c("chromosome","start","end","score")
  return(output)
}

regionwidth = 250000
nperm       = 100

realresults <-  testHetrozygousity(calldata, arrayAnnotation, regionwidth)
plot(realresults[,"score"], t = 'l')

permutations <- NULL
for(x in 1:nperm){
  arrayAnnotation[,"Probe.Set.ID"] <- arrayAnnotation[sample(nrow(arrayAnnotation)), "Probe.Set.ID"]    # Shuffle the ProbeSetIDs to new locations
  permutations <-  c(permutations, max(testHetrozygousity(calldata, arrayAnnotation, regionwidth)[,"score"]))
}

plot(permutations[sort(permutations,index.return=TRUE)$ix],type='l')
threshold <- permutations[sort(permutations,index.return=TRUE)$ix][length(permutations) * .95]


plot(realresults[,"score"], t = 'o')
abline(h=threshold)

sortperm <- permutations[sort(permutations,index.return=TRUE)$ix]
realresults <- cbind(realresults, pvalue = NA)
for(x in 1:nrow(realresults)){
  score <- realresults[x,"score"]
  realresults[x,"pvalue"] <-  1 - (sum(sortperm < realresults[x,"score"])/nperm)
}
