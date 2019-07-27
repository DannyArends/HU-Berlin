# Split in different files using:

execute <- function(x, intern = FALSE){
  cat("----", x, "\n")
  res <- system(x, intern = intern)
  cat(">>>>", res[1], "\n")
  if(res[1] >= 1) q("no")
}

for(chr in c(1:19, "X", "Y", "MT")){
  command <- paste0("samtools view -m 5G -@ 8 -b merged_sorted.bam ", chr, " > merged_chr", chr, ".bam")
  execute(command)
  command <- paste0("samtools index -@ 8 -b merged_chr", chr, ".bam merged_chr", chr, ".bai")
  execute(command)
}
q("no")