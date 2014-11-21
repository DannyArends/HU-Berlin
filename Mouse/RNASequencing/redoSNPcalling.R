fileBases <- c("5068_GAGTGG_L004_", "5069_AGTCAA_L004_",                                # BFMI samples
               "4868_GCCAAT_L001_", "5067_ATCACG_L004_",                                # B6N samples
               "5070_CGATGT_L005_", "5071_CCGTCC_L005_", "5072_TAGCTT_L005_",           # maternal B6N samples
               "5073_TTAGGC_L006_", "5074_GATCAG_L006_", "5075_ATGTCA_L006_")           # maternal BFMI samples

gatk          <- "/opt/GenomeAnalysisTK-3.2-2/GenomeAnalysisTK.jar"
referenceDir  <- "genomes"
referenceName <- "Mus_musculus.GRCm38.74.dna"
reference     <- paste0(referenceDir, "/", referenceName, ".fasta")
knownsnps     <- paste0(referenceDir, "/mgp.v3.snps.rsIDdbSNPv137.vcf")                                # Reference SNPs, download from: ftp://ftp-mouse.sanger.ac.uk/

for(fileBase in fileBases){
  outputSIRBAM  <- paste0(fileBase, "P_trimmed.aligned.sorted.realigned.recalibrated.bam")
  outputVCF     <- paste0(fileBase, ".snps.v2.vcf")
  logname       <- paste0(fileBase, ".snps.log")
  settings      <- "-stand_call_conf 20.0 -stand_emit_conf 15.0"                                         # Slightly lowered stand_emit_conf to get more DP at a SNP

  command <- paste0("nohup java -Xmx4g -jar ", gatk, " -T UnifiedGenotyper -R ", reference, " -I ", outputSIRBAM, "  --dbsnp ", knownsnps, " ", settings, " -o ", outputVCF, " > ", logname ," &")
  cat(command,"\n")
}