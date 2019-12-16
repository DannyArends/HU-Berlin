#BiocManager::install(c("MotifDb","seqLogo", "motifStack", "Biostrings", "GenomicFeatures"))

library(MotifDb)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)

TFBSjasparALL <- query(MotifDb, "Hsapiens")

mdata <- "ACACACACACACACACACACACACACACACACACACACACACACACACACACTACATAAATACAGTTTATGGGGGGTTGTTGTTTTAGAATTTGCGTGCCAGGCTGGAAGGTGGCTCAGTGTTTTTAAGAGCACTTGCTGCTCTTGCACAGGACCCAAGATTAGGTTTCTTAACAAACACATGGGAGTCCCAAGGCCCTGGCCTCCTTAGAAGGCAACACGCATGTGATGCACATTAATACATCAGGCACTATTCATATACACACAAAACAAAAATAAATATTTGAGGGGAAAATTTGTATGACAAAGCCTAAAATGAAGAACCAGGAGTAGAACAAGGCACAAAATATTCGCTTCTCCAGAGTTTCATTCAGAGCTGAAGGAACATGCATGAGCGAAGCTCCAGCGAAGCTCCAAGCAGATCGCTGGCCCAGCGGCCCAGAGCTGGGAGGCAGAGAATGCGCCCCGGTTGCCTTGGAAACCGGCCAGCGTCTGAGGACCGTTCCCTCCGCCCCCCTTTCCTGGTGCAGCTCCCCGGGAGAGGGCGGGGCTTCACAGTGAGCTAAGATGGGTCACGTGTCGAGCCCGCGTTGCTGGGCGCCGCTGCCCTGACGCAGGCCGGCTCGTCAGTCTAGTGGTGGGCTCTCCGAGGGAAGGCTTCCACAGCGGGTCCCGCCTCTCCGTGCGGCCGCCATGGATCTG"

mdata <- as.character(gsub("\n", "", mdata))

splitted <- strsplit(mdata,"")[[1]]

png("d:/tfbs.png", width = 8000, height = 1200, res = 300, pointsize = 8)

plot(c(1,length(splitted)), c(0, 1), t = 'n', xaxs='i')
for(x in 1:length(splitted)){
  text(x, 0.1, splitted[x], cex=0.4, col=as.numeric(factor(splitted[x], levels=c("A", "C", "T", "G"))))
}
abline(v = as.numeric(gregexpr("GCCAGCG", mdata)[[1]]))
abline(v = as.numeric(gregexpr("CAT", mdata)[[1]][10]))

del <- gregexpr("GCGAAGCTCCAGCGAAGCTCCA", mdata)[[1]]
abline(v = as.numeric(del), col = "red")
abline(v = as.numeric(del + .5 * attr(del, "match.length")), col = "red")
abline(v = as.numeric(del + .5 * attr(del, "match.length")), col = "green")
abline(v = as.numeric(del + attr(del, "match.length")), col = "green")


for(tfbs in 1:length(TFBSjasparALL)){
  # seqLogo(TFBSjasparALL[[tfbs]])
  TFBSjaspar <- round(100 * TFBSjasparALL[[tfbs]])
  if(ncol(TFBSjaspar) > 30) next;
  hits <- matchPWM(TFBSjaspar, mdata , "90%", with.score=TRUE)
  if(length(hits) > 0){
    scores <- mcols(hits)$score
    for(hit in 1:length(hits)){
      tfbs.id   <- names(TFBSjasparALL)[tfbs]
      tfbs.name   <- as.character(as.data.frame(values(TFBSjasparALL [tfbs]))["geneSymbol"])
      tfbs.seq    <- as.character(unlist(hits[hit]))
      tfbs.start  <- hits@ranges@start[hit]
      tfbs.width  <- hits@ranges@width[hit]
      seqsplit <- rownames(TFBSjaspar)[apply(TFBSjaspar,2,which.max)]#strsplit(tfbs.seq, "")[[1]]
      ch <- 1
      for(l in tfbs.start:(tfbs.start+length(seqsplit) - 1)){
        text(l, 0.2 + (tfbs %% 16) / 20, seqsplit[ch], cex=0.4, col= (1+ as.numeric(seqsplit[ch] != splitted[l])))
        ch <- ch + 1
      }
      #text(tfbs.start, 0.2 + (tfbs %% 16) / 20, tfbs.seq, cex=0.4)
      cat(tfbs, tfbs.id, scores[hit], tfbs.name, tfbs.seq, tfbs.start, tfbs.width, "\n")
    }
  }
}
dev.off()


