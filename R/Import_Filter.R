#' Imports SNP data from GATK VariablesToTable output
#'
#' Imports SNP data from the output of the
#' \href{https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php}{VariantsToTable}
#' function in GATK. After importing the data, the function then calculates
#' total reference allele frequency for both bulks together, the delta SNP index
#' (i.e. SNP index of the low bulk substacted from the SNP index of the high
#' bulk), the G statistic and returns a data frame. The required GATK fields
#' (-F) are CHROM (Chromosome) and POS (Position). The required Genotype fields
#' (-GF) are AD (Allele Depth), DP (Depth), GQ  (Genotype Quality). Recommended
#' fields are REF (Reference allele) and ALT (Alternative allele) Recommended
#' Genotype feilds are PL (Phred-scaled likelihoods)
#'
#' @param filename The name of the GATK VariablesToTable output .table file which the
#'   data are to be read from.
#' @param highBulk The sample name of the High Bulk
#' @param lowBulk The sample name of the Low Bulk
#' @param chromList a string vector of the chromosomes to be used in the
#'   analysis. Useful for filtering out unwanted contigs etc.
#' @return Returns a data frame containing columns for Read depth (DP),
#'   Reference Allele Depth (AD.REF) and Alternative Allele Depth (AD.ALT),
#'   Genoytype Quality (GQ) and SNPindex for each bulk (indicated by .HIGH and
#'   .LOW column name suffix). Total reference allele frequnce "REF_FRQ" is the
#'   sum of AD.REF for both bulks divided by total Depth for that SNP. The
#'   deltaSNPindex is equal to  SNPindex.HIGH - SNPindex.LOW. The GStat column
#'   is the calculated G statistic for that SNP.
#' @seealso \code{\link{GetGStat}} for explaination of how G statistic is
#'   calculated.
#'   \href{http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it}{What
#'   is a VCF and how should I interpret it?} for more information on GATK
#'   Fields and Genotype Fields
#' @examples df <-  ImportFromGATK(filename = file.table,
#'     HighBulk = HighBulkSampleName,
#'     LowBulk = LowBulkSampleName,
#'     ChromList = c("Chr1","Chr4","Chr7"))
#' @export importFromGATK

importFromGATK <- function(filename,
    highBulk = character(),
    lowBulk = character(),
    chromList = NULL) {
    message("Importing SNPs from file")
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
    
    #Keep only wanted chromosomes
    if (!is.null(ChromList)) {
        SNPset <- SNPset[SNPset$CHROM %in% ChromList, ]
    }
    #arrange the chromosomes by natural order sort, eg Chr1, Chr10, Chr2 >>> Chr1, Chr2, Chr10
    SNPset$CHROM <- factor(SNPset$CHROM, levels = gtools::mixedsort(unique(SNPset$CHROM)))
    
    # Calculate some descriptors
    SNPset$REF_FRQ <-
        (SNPset$AD_REF.HIGH + SNPset$AD_REF.LOW) / (SNPset$DP.HIGH + SNPset$DP.LOW)
    SNPset$deltaSNP <- SNPset$SNPindex.HIGH - SNPset$SNPindex.LOW
    return(SNPset)
}


#' Filter SNPs based on read depth and quality
#'
#' Use filtering paramaters to filter out high and low depth reads as well as
#' low Genotype Quality as defined by GATK. All filters are optional but recommended.
#'
#' @param SNPset The data frame imported by \code{ImportFromGATK}
#' @param refAlleleFreq A numeric < 1. This will filter out SNPs with a
#'   Reference Allele Frequency less than \code{refAlleleFreq} and greater than
#'   1 - \code{refAlleleFreq}. Eg. \code{refAlleleFreq = 0.3} will keep SNPs
#'   with 0.3 <= REF_FRQ <= 0.7
#' @param filterAroundMedianDepth Filters total SNP read depth for both bulks. A
#'   median and median absolute deviation (MAD) of depth will be calculated.
#'   SNPs with read depth greater or less than \code{filterAroundMedianDepth}
#'   MADs away from the median will be filtered.
#' @param minTotalDepth The minimum total read depth for a SNP (counting both
#'   bulks)
#' @param maxTotalDepth The maximum total read depth for a SNP (counting both
#'   bulks)
#' @param minSampleDepth The minimum read depth for a SNP in each bulk
#' @param minGQ The minimum Genotype Quality as set by GATK. This is a measure
#'   of how confident GATK was with the assigned genotype (i.e. homozygous ref,
#'   heterozygous, homozygous alt). See
#'   \href{http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it}{What
#'   is a VCF and how should I interpret it?}
#' @param verbose logical. If \code{TRUE} will report number of SNPs filtered in
#'   each step.
#' @return Returns a subset of the data frame supplied which meets the filtering
#'   conditions applied by the selected parameters. If \code{verbose} is
#'   \code{TRUE} the function reports the number of SNPs filtered in each step
#'   as well as the initiatl number of SNPs, the total number of SNPs filtered
#'   and the remaining number.
#'
#' @seealso See \code{\link[stats]{mad}} for explaination of calculation of
#'   median absolute deviation.
#'   \href{http://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it}{What
#'   is a VCF and how should I interpret it?} for more information on GATK
#'   Fields and Genotype Fields
#' @examples df_filt <- FilterSNPs(
#'     df,
#'     refAlleleFreq = 0.3,
#'     minTotalDepth = 40,
#'     maxTotalDepth = 80,
#'     minSampleDepth = 20,
#'     minGQ = 99,
#'     verbose = TRUE
#' )
#'
#' @export filterSNPs

filterSNPs <- function(SNPset,
    refAlleleFreq = NULL,
    filterAroundMedianDepth = 2.5,
    minTotalDepth,
    maxTotalDepth,
    minSampleDepth = NULL,
    minGQ = 99,
    verbose = TRUE) {
    org_count <- nrow(SNPset)
    count <- nrow(SNPset)
    
    # Filter by total reference allele frequency
    if (!is.null(refAlleleFreq)) {
        if (verbose) {
            message(
                "Filtering by reference allele frequency: ",
                refAlleleFreq,
                " <= REF_FRQ <= ",
                1 - refAlleleFreq
            )
        }
        SNPset <- dplyr::filter(SNPset, SNPset$REF_FRQ < 1 - refAlleleFreq &
                SNPset$REF_FRQ > refAlleleFreq)
            if (verbose) {
            message("...Filtered ", count - nrow(SNPset), " SNPs")
        }
        count <- nrow(SNPset)
    }
    
    #Total read depth filtering
    
    if (!missing(filterAroundMedianDepth)) {
        # filter by Read depth for each SNP FilterByMAD MADs around the median
        madDP <-
            mad(
                x = (SNPset$DP.HIGH + SNPset$DP.LOW),
                constant = 1,
                na.rm = TRUE
            )
        medianDP <-
            median(x = (SNPset$DP.HIGH + SNPset$DP.LOW),
                na.rm = TRUE)
        maxDP <- medianDP + filterAroundMedianDepth * madDP
        minDP <- medianDP - filterAroundMedianDepth * madDP
        SNPset <- dplyr::filter(SNPset, (DP.HIGH + DP.LOW) <= maxDP &
                (DP.HIGH + DP.LOW) >= minDP)

        if(verbose) {message("Filtering by total read depth: ",
            filterAroundMedianDepth,
            " MADs arround the median: ", minDP, " <= Total DP <= ", maxDP)
            message("...Filtered ", count - nrow(SNPset), " SNPs")}
        count <- nrow(SNPset)
        
    }
    
    if (!missing(minTotalDepth)) {
        # Filter by minimum total SNP depth
        if (verbose) {
            message("Filtering by total sample read depth: Total DP >= ",
                minTotalDepth)
        }
        SNPset <-
            dplyr::filter(SNPset, (DP.HIGH + DP.LOW) >= minTotalDepth)

        if (verbose) {
            message("...Filtered ", count - nrow(SNPset), " SNPs")
        }
        count <- nrow(SNPset)
    }
    
    if (!missing(maxTotalDepth)) {
        # Filter by maximum total SNP depth
        if (verbose) {
            message("Filtering by total sample read depth: Total DP <= ",
                maxTotalDepth)
        }
        SNPset <-
            dplyr::filter(SNPset, (DP.HIGH + DP.LOW) <= maxTotalDepth)
        if (verbose) {
            message("...Filtered ", count - nrow(SNPset), " SNPs")
        }
        count <- nrow(SNPset)
    }
    
    
    # Read depth in each bulk should be greater than 40
    if (!missing(minSampleDepth)) {
        if (verbose) {
            message("Filtering by per sample read depth: DP >= ",
                minSampleDepth)
        }
        SNPset <-
            dplyr::filter(SNPset, DP.HIGH >= minSampleDepth &
                    SNPset$DP.LOW >= minSampleDepth)
        if (verbose) {
            message("...Filtered ", count - nrow(SNPset), " SNPs")
        }
        count <- nrow(SNPset)
    }
    
    # Filter by LOW BULK Genotype Quality
    if (!is.null(minGQ)) {
        if (verbose) {
            message("Filtering by Genotype Quality: GQ >= ", minGQ)
        }
        SNPset <-
            dplyr::filter(SNPset, GQ.LOW >= minGQ & GQ.HIGH >= minGQ)
        if (verbose) {
            message("...Filtered ", count - nrow(SNPset), " SNPs")
        }
        count <- nrow(SNPset)
    }
    
    # #Filter SNP Clusters
    # if (!is.null(SNPsInCluster) & !is.null(ClusterWin)) {
    #     tmp <- which(diff(SNPset$POS, SNPsInCluster-1) < ClusterWin)
    # message("...Filtered ", count - nrow(SNPset), " SNPs")
    # count <- nrow(SNPset)
    # }
    if (verbose) {
        message(
            "Original SNP number: ",
            org_count,
            ", Filtered: ",
            org_count - count,
            ", Remaining: ",
            count
        )
    }
    return(SNPset)
}
