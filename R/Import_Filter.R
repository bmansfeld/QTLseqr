# Imports SNP data from GATK VarToTable output.
# The required GATK fields (-F) are CHROM (Chromosome) and POS (Position)
# The required Genotype fields (-GF) are AD (Allele Depth), DP (Depth), GQ  (Genotype Quality)
# Recommended fields are REF (Reference allele) and ALT (Alternative allele)
# Recommended Genotype feilds are PL (Phred-scaled likelihoods)
# After importing the data, the function then calculates total reference allele frequency for both bulks together,
# the delta SNP index (i.e. SNP index of the low bulk substracted from the SNP index of the high bulk)
# and the G statistic

ImportFromGATK <- function(filename,
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
    if (!is.null(ChromList)) {
        SNPset <- subset(SNPset, CHROM %in% ChromList)
    }
    # Calculate some descriptors
    SNPset$REF_FRQ <-
        (SNPset$AD_REF.HIGH + SNPset$AD_REF.LOW) / (SNPset$DP.HIGH + SNPset$DP.LOW)
    SNPset$deltaSNP <- SNPset$SNPindex.HIGH - SNPset$SNPindex.LOW

    # calculate G Statistic
    SNPset$GStat <- GetGStat(SNPset)
    return(SNPset)
}


# Filter SNPs based on some usefull parameters including read depth and quality
FilterSNPs <- function(SNPset,
    RefAlleleFreq = NULL,
    FilterByMAD = 2.5,
    MinTotalDepth,
    MaxTotalDepth,
    MinSampleDepth = NULL,
    MinGQ = 99) {
    # Filter by total reference allele frequency
    if (!is.null(RefAlleleFreq)) {
        message(
            "Filtering by reference allele frequency: Min = ",
            RefAlleleFreq,
            ", Max = ",
            1 - RefAlleleFreq
        )
        SNPset <-
            subset(SNPset,
                REF_FRQ < 1 - RefAlleleFreq & REF_FRQ > RefAlleleFreq)
    }

    #Total read depth filtering

    if (!is.null(FilterByMAD)) {
        # filter by Read depth for each SNP FilterByMAD MADs around the median
        # constant is set at 1.4826 assuming reads have normal distribution
        message("Filtering by total read depth: ",
            FilterByMAD,
            " MADs arround the median")
        madDP <-
            mad(x = (SNPset$DP.HIGH + SNPset$DP.LOW),
                constant = 1.4826)
        medianDP <- median(SNPset$DP.HIGH + SNPset$DP.LOW)
        SNPset <-
            subset(
                SNPset,
                (DP.HIGH + DP.LOW) < (medianDP + FilterByMAD * madDP) &
                    (DP.HIGH + DP.LOW) > (medianDP - FilterByMAD * madDP)
            )

    }

    if (!missing(MinTotalDepth)) {
        # Filter by minimum total SNP depth
        SNPset <- subset(SNPset, (DP.HIGH + DP.LOW) > MinTotalDepth)
    }

    if (!missing(MaxTotalDepth)) {
        # Filter by maximum total SNP depth
        SNPset <- subset(SNPset, (DP.HIGH + DP.LOW) < MaxTotalDepth)
    }


    # Read depth in each bulk should be greater than 40
    if (!is.null(MinSampleDepth)) {
        SNPset <-
            subset(SNPset,
                DP.HIGH >= MinSampleDepth & DP.LOW >= MinSampleDepth)
    }

    # Filter by LOW BULK Genotype Quality
    if (!is.null(MinGQ)) {
        SNPset <-
            SNPset <-
            subset(SNPset, GQ.LOW == MinGQ & GQ.HIGH == MinGQ)
    }

    return(SNPset)
}
