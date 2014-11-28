
growth <- phenotypes[F2, c("d21", "d28", "d35", "d42", "d49", "d56", "d63", "d70", "d71")]

plot(c(1, 9), c(0,max(growth,na.rm=TRUE)), t='n', ylab="Bodyweight")
apply(growth,1,function(x){ points(x, t='l') })

aa <- which(abs(apply(growth[,c(8:9)],1,function(x){x[2] - x[1]})) > 3)       # Which individuals grew more then 3 grams in the last day ?
# 6661109 6661116 6661339 6661340 6661341 6661342 6661343 6662008

mri <- phenotypes[F2, c("mri42d_fat", "mri42d_lean", "mri56d_fat", "mri56d_lean", "mri70d_fat", "mri70d_lean")]
mri[which(mri[,5] < 0), ]                                               # Which individuals have negative fat weight ?
# 6661341

# Fat mass
plot(c(1, 3), c(0,max(mri[,c(1, 3, 5)], na.rm=TRUE)), t='n', ylab="FAT")
apply(mri[,c(1, 3, 5)],1,function(x){ points(x, t='l', col=2+ sign(cor(1:3, as.numeric(x)))) })
updown <- apply(mri[,c(1, 3, 5)],1,function(x){ sign(cor(1:3, as.numeric(x))) })
#661105 6661107 6661108 6661124 6661156 6661227 6661261 6661308 6661326 6661340 6661341 6661342 6661343 6661386 6661389 6661408 6661448 6661456 6661457 6661460 6661554 6661561 6661758 6661817 
#6662171 6662172 6662173 6662243 6662246 


# Lean mass
plot(c(1, 3), c(0,max(mri[,c(2, 4, 6)], na.rm=TRUE)), t='n', ylab="LEAN")
apply(mri[,c(2, 4, 6)],1,function(x){ points(x, t='l', col=2+ sign(cor(1:3, as.numeric(x)))) })
updown <- apply(mri[,c(2, 4, 6)],1,function(x){ sign(cor(1:3, as.numeric(x))) })
# 6661339 6661341 

