# Find a minimal SET of SNPs describing the popualtion
  
table(unlist(absnpdata["snp993-scaffold1026-156620",rownames(samples)[which(samples[,"Breed"] == "Nu")]]))
table(unlist(absnpdata["snp993-scaffold1026-156620",rownames(samples)[which(samples[,"Breed"] == "Ni")]]))
table(unlist(absnpdata["snp993-scaffold1026-156620",rownames(samples)[which(samples[,"Breed"] == "Dese")]]))
table(unlist(absnpdata["snp993-scaffold1026-156620",rownames(samples)[which(samples[,"Breed"] == "Tagg")]]))


table(unlist(absnpdata["snp36346-scaffold4353-5630",rownames(samples)[which(samples[,"Breed"] == "Nu")]]))
table(unlist(absnpdata["snp36346-scaffold4353-5630",rownames(samples)[which(samples[,"Breed"] == "Ni")]]))
table(unlist(absnpdata["snp36346-scaffold4353-5630",rownames(samples)[which(samples[,"Breed"] == "Dese")]]))
table(unlist(absnpdata["snp36346-scaffold4353-5630",rownames(samples)[which(samples[,"Breed"] == "Tagg")]]))

Nu.ind <- rownames(samples)[which(samples[,"Breed"] == "Nu")]
Tagg.ind <- rownames(samples)[which(samples[,"Breed"] == "Tagg")]

Nu.specs <- which(lapply(apply(numsnpdata[,rownames(samples)[which(samples[,"Breed"] == "Nu")]],1,table), length) == 1)
Nu.geno <- numsnpdata[Nu.specs,]

Ta.specs <- which(lapply(apply(numsnpdata[,rownames(samples)[which(samples[,"Breed"] == "Tagg")]],1,table), length) == 1)
Ta.geno <- numsnpdata[Ta.specs,]


  Nu.geno[,Nu.ind][1:5,1:5]
  Nu.geno[,Tagg.ind][1:5,1:5]
  Ta.geno[,Nu.ind][1:5,1:5]
  Ta.geno[,Tagg.ind][1:5,1:5]

Nu.geno["snp2670-scaffold1077-230230",Nu.ind]
Nu.geno["snp2670-scaffold1077-230230",Tagg.ind]

length(which(
  Nu.geno["snp2670-scaffold1077-230230",Tagg.ind] == 3 & 
  Nu.geno["snp55670-scaffold863-1403541",Tagg.ind] == 1 &
  Nu.geno["snp37429-scaffold456-4388525",Tagg.ind] == 3 &
  Nu.geno["snp58374-scaffold947-3741069",Tagg.ind] == 1 &
  Nu.geno["snp58287-scaffold947-140466",Tagg.ind] == 3 &
  Nu.geno["snp31510-scaffold349-2559367",Tagg.ind] == 1 &
  Nu.geno["snp57323-scaffold912-2798776",Tagg.ind] == 1 &
  Nu.geno["snp33653-scaffold395-1378592",Tagg.ind] == 1
)) / length(Tagg.ind)


length(which(Nu.geno["snp2670-scaffold1077-230230",Tagg.ind] == 3 & Nu.geno["snp55670-scaffold863-1403541",Tagg.ind] == 1)) / length(Tagg.ind)
length(which(Nu.geno["snp2670-scaffold1077-230230",Tagg.ind] == 3 & Nu.geno["snp55670-scaffold863-1403541",Tagg.ind] == 1)) / length(Tagg.ind)


