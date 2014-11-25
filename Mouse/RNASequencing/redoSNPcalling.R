fileBases <- c("5068_GAGTGG_L004_", "5069_AGTCAA_L004_",                                # BFMI samples
               "4868_GCCAAT_L001_", "5067_ATCACG_L004_",                                # B6N samples
               "5070_CGATGT_L005_", "5071_CCGTCC_L005_", "5072_TAGCTT_L005_",           # maternal B6N samples
               "5073_TTAGGC_L006_", "5074_GATCAG_L006_", "5075_ATGTCA_L006_")           # maternal BFMI samples

gatk          <- "/opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"
referenceDir  <- "genomes"
referenceName <- "Mus_musculus.GRCm38.74.dna"
reference     <- paste0(referenceDir, "/", referenceName, ".fasta")
knownsnps     <- paste0(referenceDir, "/mgp.v3.snps.rsIDdbSNPv137.vcf")                                # Reference SNPs, download from: ftp://ftp-mouse.sanger.ac.uk/
mpileup       <- "~/Github/samtools/samtools mpileup"

for(fileBase in fileBases){
  outputSIRBAM  <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam")
  outputBCF     <- paste0(fileBase, ".bcf")
  outputVarBCF  <- paste0(fileBase, "var.bcf")
  outputVCF     <- paste0(fileBase, ".snps.bcftools.vcf")
  logname       <- paste0(fileBase, ".snps.log")

  command <- paste0("nohup ", mpileup, " -g -f ",reference," ", outputSIRBAM, " -o ", outputBCF, " &")
  cat(command,"\n")
  command <- paste0("nohup bcftools call -c ", outputBCF, " | ~/Github/bcftools/vcfutils.pl varFilter -d 10 - > ", outputVCF, " &")
  cat(command,"\n")
}
