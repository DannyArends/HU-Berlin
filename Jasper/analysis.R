#
# Jasper analysis pig data
#
#


setwd("d:/Jasper")

## CD 3 - 7mm

mdata <- read.table("CD3.data.txt",sep="\t",header=TRUE, row.names=1)

for(x in 2:7){
  pv1 <- t.test(mdata[mdata[,"Group"] == 1, x], mdata[mdata[,"Group"] == 2, x])
  pv2 <- t.test(mdata[mdata[,"Group"] == 2, x], mdata[mdata[,"Group"] == 3, x])
  pv3 <- t.test(mdata[mdata[,"Group"] == 1, x], mdata[mdata[,"Group"] == 3, x])
  cat(colnames(mdata)[x], ":", pv1$p.value, pv2$p.value, pv3$p.value, "\n")
}


mdata <- read.table("CD3.comparison.txt",sep="\t",header=TRUE)
anova(lm(mdata[,"Gesamt"] ~ mdata[,"Group"] + mdata[,"location"]))
anova(lm(mdata[,"CD3"] ~ mdata[,"Group"] + mdata[,"location"]))

# Conclusions:
Combined Rand and Zentrum
  - Gesamt measurement: Significant difference between locations: Rand versus Zentrum (P = 0.02398)
  - Gesamt measurement : Significant difference between groups: 1, 2 and 3 (P = 0.08988)

  - CD3 measurement:  Highly Significant difference between locations: Rand versus Zentrum (P = 2.833e-10)
  - CD3 measurement : Significant difference between groups: 1, 2 and 3 (P = 0.0608)

- Single groups:

Keimzentrum.gesamt :      0.06082025.  0.534938  0.20149 
Keimzentrum.CD3 :         0.1557976   0.6617248 0.1401711 
Rand.gesamt :             0.2645785   0.9686234 0.3176751 
Rand.CD3 :                0.02705842* 0.2712357 0.03868401*
Paracortex.gesamt :       0.7522009   0.8958058 0.8237213 
Paracortes.CD3 :          0.6903786   0.5806389 0.9534926 


## CD 4 - 7mm

mdata <- read.table("CD4.data.txt",sep="\t",header=TRUE, row.names=1)

for(x in 2:7){
  pv1 <- t.test(mdata[mdata[,"Group"] == 1, x], mdata[mdata[,"Group"] == 2, x])
  pv2 <- t.test(mdata[mdata[,"Group"] == 2, x], mdata[mdata[,"Group"] == 3, x])
  pv3 <- t.test(mdata[mdata[,"Group"] == 1, x], mdata[mdata[,"Group"] == 3, x])
  cat(colnames(mdata)[x], ":", pv1$p.value, pv2$p.value, pv3$p.value, "\n")
}
mdata <- read.table("CD4.comparison.txt",sep="\t",header=TRUE)
anova(lm(mdata[,"Gesamt"] ~ mdata[,"Group"] + mdata[,"location"]))
anova(lm(mdata[,"CD4"] ~ mdata[,"Group"] + mdata[,"location"]))

# Conclusions:
Combined Rand and Zentrum
  - Gesamt measurement: NO Significant difference between locations: Rand versus Zentrum (P = 0.74364)
  - Gesamt measurement : Significant difference between groups: 1, 2 and 3 (P = 0.06467)
  
  - CD4 measurement:  Highly Significant difference between locations: Rand versus Zentrum (P = 0.000192)
  - CD4 measurement : No Significant difference between groups: 1, 2 and 3 (P = 0.286769)
  
                              1-2       2-3       1-3
    Keimzentrum.gesamt :    0.1591979 0.1931342 0.3441863 
    Keimzentrum.CD4 :       0.4293258 0.9608805 0.4443296 
    Rand.gesamt :           0.1846171 0.9146737 0.194292 
    Rand.CD4 :              0.6478391 0.7207282 0.8291717 
    Paracortex.gesamt :     0.4211689 0.4983198 0.6286881 
    Paracortes.CD4 :        0.6636186 0.7143741 0.8233504 


## CD 8a - 7mm

mdata <- read.table("CD8a.data.txt",sep="\t",header=TRUE, row.names=1)

for(x in 2:7){
  pv1 <- t.test(mdata[mdata[,"Group"] == 1, x], mdata[mdata[,"Group"] == 2, x])
  pv2 <- t.test(mdata[mdata[,"Group"] == 2, x], mdata[mdata[,"Group"] == 3, x])
  pv3 <- t.test(mdata[mdata[,"Group"] == 1, x], mdata[mdata[,"Group"] == 3, x])
  cat(colnames(mdata)[x], ":", pv1$p.value, pv2$p.value, pv3$p.value, "\n")
}

mdata <- read.table("CD8a.comparison.txt",sep="\t",header=TRUE)
anova(lm(mdata[,"Gesamt"] ~ mdata[,"Group"] + mdata[,"location"]))
anova(lm(mdata[,"CD8a"] ~ mdata[,"Group"] + mdata[,"location"]))


# Conclusions:
Combined Rand and Zentrum
  - Gesamt measurement: No Significant difference between locations: Rand versus Zentrum (P = 0.9682)
  - Gesamt measurement : No Significant difference between groups: 1, 2 and 3 (P = 0.1221)

  - CD8a measurement:  Highly Significant difference between locations: Rand versus Zentrum (P = 0.0003082)
  - CD8a measurement : No Significant difference between groups: 1, 2 and 3 (P = 0.1894574)
  
                              1-2       2-3       1-3
    Keimzentrum.gesamt :  0.0713338.   0.4603786  0.1981134 
    Keimzentrum.CD8a :    0.05088885.  0.2596316  0.2802181 
    Rand.gesamt :         0.02281575*  0.3997435  0.1953908 
    Rand.CD8a :           0.5197883    0.696515   0.4267687 
    Paracortex.gesamt :   0.4083458    0.4029569  0.8934307 
    Paracortes.CD8a :     0.72508      0.4487628  0.3236331 
    
    