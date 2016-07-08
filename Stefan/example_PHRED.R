##### replace phredscore with
### get special characters
phred <-  unlist(strsplit(rawToChar(as.raw(33:126)),""))
phrednumeric <- c(1:94)
aaa <- cbind(phred,phrednumeric)
aaa <- aaa[1:42,]


seqs <-
"AAA1AA1AF3D3BAA3EEEFFGGHFCBFHHHFGHCFFF1FEGG2FFGHAE1BFEE0FGHGFBAAGH1FE//EFF010F0AG/FFH2BF1FECHFGHFEFHEGG?FF1BDFFA0B2BE/>FC2B1BGGBFG<1<ECCC//?FF0?1F1?11"

seqs <- strsplit(seqs, "")

# replace seqs w/ rownumbers of aaa
transformed <- unlist(lapply(unlist(seqs), function(x){
           return(aaa[match(x,aaa[,1]),2])
}))

names(transformed) <- unlist(seqs)

