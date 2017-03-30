
# Imports SNP data from GATK VarToTable output.
# The required GATK fields (-F) are CHROM (Chromosome) and POS (Position)
# The required Genotype fields (-GF) are AD (Allele Depth), DP (Depth), GQ  (Genotype Quality)
# Recommended fields are REF (Reference allele) and ALT (Alternative allele)
# Recommended Genotype feilds are PL (Phred-scaled likelihoods)
# After importing the data, the function then calculates total reference allele frequency for both bulks together,
# the delta SNP index (i.e. SNP index of the low bulk substracted from the SNP index of the high bulk)
# and the G statistic
import_GATK_VCF <- function(filename,
                            HighBulk = character(),
                            LowBulk = character(),
                            ChromList = NULL) {
  VarTable <-
    read.table(file = filename,
               header = T,
               stringsAsFactors = F)

   # Format data frame for analysis
  SNPset <- VarTable[, 1:4]

  # High Bulk data
  SNPset$DP.HIGH <- VarTable[, paste0(HighBulk, ".DP")]
  SNPset$AD_REF.HIGH <-
    as.numeric(gsub(",.*$", "", x = VarTable[, paste0(HighBulk, ".AD")]))
  SNPset$AD_ALT.HIGH <- SNPset$DP.HIGH - SNPset$AD_REF.HIGH
  SNPset$GQ.HIGH <- VarTable[, paste0(HighBulk, ".GQ")]
  # Calculate SNP index
  SNPset$SNPindex.HIGH <- SNPset$AD_ALT.HIGH / SNPset$DP.HIGH

  # Low Bulk data
  SNPset$DP.LOW <- VarTable[, paste0(LowBulk, ".DP")]
  SNPset$AD_REF.LOW <-
    as.numeric(gsub(",.*$", "", x = VarTable[, paste0(LowBulk, ".AD")]))
  SNPset$AD_ALT.LOW <- SNPset$DP.LOW - SNPset$AD_REF.LOW
  SNPset$GQ.LOW <- VarTable[, paste0(LowBulk, ".GQ")]
  SNPset$SNPindex.LOW <- SNPset$AD_ALT.LOW / SNPset$DP.LOW

  #Subset any unwanted chromosomes
  SNPset <- subset(SNPset, CHROM %in% ChromList)

  # Calculate some descriptors
  SNPset$REF_FRQ <- (SNPset$AD_REF.HIGH + SNPset$AD_REF.LOW) / (SNPset$DP.HIGH + SNPset$DP.LOW)
  SNPset$deltaSNP <- SNPset$SNPindex.HIGH - SNPset$SNPindex.LOW

  # calculate G Statistic
  SNPset$GStat <- GetGStat(SNPset)
  return(SNPset)
}
