#import vcf from GATK VarToTable
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

  # Calculate some descripters
  SNPset$REF_FRQ <- (SNPset$AD_REF.HIGH + SNPset$AD_REF.LOW) / (SNPset$DP.HIGH + SNPset$DP.LOW)
  SNPset$deltaSNP <- SNPset$SNPindex.HIGH - SNPset$SNPindex.LOW
  return(SNPset)
}
