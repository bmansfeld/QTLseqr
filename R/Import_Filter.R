#' Imports SNP data from GATK VariablesToTable output
#'
#' Imports SNP data from the output of the
#' \href{https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php}{VariantsToTable}
#' function in GATK. After importing the data, the function then calculates
#' total reference allele frequency for both bulks together, the delta SNP index
#' (i.e. SNP index of the low bulk substacted from the SNP index of the high
#' bulk), the G statistic and returns a data frame. The required GATK fields
#' (-F) are CHROM (Chromosome) and POS (Position). The required Genotype fields
#' (-GF) are AD (Allele Depth), DP (Depth). Recommended
#' fields are REF (Reference allele) and ALT (Alternative allele) Recommended
#' Genotype feilds are PL (Phred-scaled likelihoods) and GQ  (Genotype Quality).
#'
#' @param file The name of the GATK VariablesToTable output .table file which the
#'   data are to be read from.
#' @param highBulk The sample name of the High Bulk
#' @param lowBulk The sample name of the Low Bulk
#' @param chromList a string vector of the chromosomes to be used in the
#'   analysis. Useful for filtering out unwanted contigs etc.
#' @return Returns a data frame containing columns for Read depth (DP),
#'   Reference Allele Depth (AD_REF) and Alternative Allele Depth (AD_ALT),
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
#'     highBulk = highBulkSampleName,
#'     lowBulk = lowBulkSampleName,
#'     chromList = c("Chr1","Chr4","Chr7"))
#' @export importFromGATK

importFromGATK <- function(file,
    highBulk = character(),
    lowBulk = character(),
    chromList = NULL) {
    
    SNPset <-
        readr::read_tsv(file = file, 
                        col_names = TRUE, 
                        col_types = readr::cols(.default = readr::col_guess(),
                                         CHROM = "c", POS = "i")
                        )
    
    if (!all(
        c(
        "CHROM", 
        "POS", 
        paste0(highBulk, ".AD"), 
        paste0(lowBulk, ".AD"), 
        paste0(highBulk, ".DP"), 
        paste0(lowBulk, ".DP")
        ) %in% names(SNPset))) {
        stop("One of the required fields is missing. Check your table file.")
    }
    
# rename columns based on bulk names and flip headers (ie HIGH.AD -> AD.HIGH to match the rest of the functions
    colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")] <- 
        gsub(pattern = highBulk, 
             replacement = "HIGH", 
             x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")])
    colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")] <- 
        gsub(pattern = lowBulk, 
             replacement = "LOW", 
             x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")])
    
    colnames(SNPset) <-
        sapply(strsplit(colnames(SNPset), "[.]"),
            function(x) {paste0(rev(x),collapse = '.')})
    
    #Keep only wanted chromosomes
    if (!is.null(chromList)) {
        message("Removing the following chromosomes: ", paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in% chromList], collapse = ", "))
        SNPset <- SNPset[SNPset$CHROM %in% chromList, ]
    }
    #arrange the chromosomes by natural order sort, eg Chr1, Chr10, Chr2 >>> Chr1, Chr2, Chr10
    SNPset$CHROM <-
        factor(SNPset$CHROM, levels = gtools::mixedsort(unique(SNPset$CHROM)))
    
    SNPset <- SNPset %>%
        tidyr::separate(
            col = AD.LOW,
            into = "AD_REF.LOW",
            sep = ",",
            extra = "drop",
            convert = TRUE
        ) %>%
        tidyr::separate(
            col = AD.HIGH,
            into = "AD_REF.HIGH",
            sep = ",",
            extra = "drop",
            convert = TRUE
        ) %>%
        dplyr::mutate(
            AD_ALT.HIGH = DP.HIGH - AD_REF.HIGH,
            AD_ALT.LOW = DP.LOW - AD_REF.LOW,
            SNPindex.HIGH = AD_ALT.HIGH / DP.HIGH,
            SNPindex.LOW = AD_ALT.LOW / DP.LOW,
            REF_FRQ = (AD_REF.HIGH + AD_REF.LOW) / (DP.HIGH + DP.LOW),
            deltaSNP = SNPindex.HIGH - SNPindex.LOW
        ) %>%
        dplyr::select(
            -dplyr::contains("HIGH"),
            -dplyr::contains("LOW"),
            -dplyr::one_of("deltaSNP", "REF_FRQ"),
            dplyr::matches("AD.*.LOW"),
            dplyr::contains("LOW"),
            dplyr::matches("AD.*.HIGH"),
            dplyr::contains("HIGH"),
            dplyr::everything()
        )
    
    return(as.data.frame(SNPset))
}


#' Import SNP data from a delimited file
#'
#' After importing the data from a delimited file, the function then calculates
#' total reference allele frequency for both bulks together, the delta SNP index
#' (i.e. SNP index of the low bulk substacted from the SNP index of the high
#' bulk), the G statistic and returns a data frame. The required columns in 
#' the file are CHROM (Chromosome) and POS (Position) as well as the refernce and alternate 
#' allele depths (number of reads supporting each allele). The allele depths should be in columns 
#' named in this format: \code{AD_(<ALT/REF>).<sampleName>}. For example, the column for alternate 
#' allele depth for a high bulk sample named "sample1", should be "AD_ALT.sample1". 
#'
#' @param file The name of the file which the
#'   data are to be read from.
#' @param highBulk The sample name of the High Bulk. Defaults to "HIGH"
#' @param lowBulk The sample name of the Low Bulk. Defaults to "LOW"
#' @param chromList a string vector of the chromosomes to be used in the
#'   analysis. Useful for filtering out unwanted contigs etc.
#' @param sep the field separator character. Values on each line of the file are separated by this character. Default is for csv file ie ",".
#' @return Returns a data frame containing columns for Read depth (DP),
#'   Reference Allele Depth (AD_REF) and Alternative Allele Depth (AD_ALT),
#'   any other SNP associated columns in the file, and SNPindex for each bulk (indicated by .HIGH and
#'   .LOW column name suffix). Total reference allele frequnce "REF_FRQ" is the
#'   sum of AD_REF for both bulks divided by total Depth for that SNP. The
#'   deltaSNPindex is equal to  SNPindex.HIGH - SNPindex.LOW.
#'
#' @export importFromTable

importFromTable <-
    function(file,
             highBulk = "HIGH",
             lowBulk = "LOW",
             chromList = NULL,
             sep = ",") {
        SNPset <-
            readr::read_delim(
                file = file,
                delim = sep,
                col_names = TRUE,
                col_types = readr::cols(
                    .default = readr::col_guess(),
                    CHROM = "c",
                    POS = "i"
                )
            )
        # check CHROM
        if (!"CHROM" %in% names(SNPset)) {
            stop("No 'CHROM' coloumn found.")
        }
        
        # check POS
        if (!"POS" %in% names(SNPset)) {
            stop("No 'POS' coloumn found.")
        }
        
        # check AD_REF.HIGH
        if (!paste0("AD_REF.", highBulk) %in% names(SNPset)) {
            stop(
                "No High Bulk AD_REF coloumn found. Column should be named 'AD_REF.highBulkName'."
            )
        }
        
        # check AD_REF.LOW
        if (!paste0("AD_REF.", lowBulk) %in% names(SNPset)) {
            stop("No Low Bulk AD_REF coloumn found. Column should be named 'AD_REF.lowBulkName'.")
        }
        
        #check AD_ALT.HIGH
        if (!paste0("AD_ALT.", highBulk) %in% names(SNPset)) {
            stop(
                "No High Bulk AD_REF coloumn found. Column should be named 'AD_ALT.highBulkName'."
            )
        }
        
        # check AD_ALT.LOW
        if (!paste0("AD_ALT.", lowBulk) %in% names(SNPset)) {
            stop("No Low Bulk AD_ALT coloumn found. Column should be named 'AD_ALT.lowBulkName'.")
        }
        
        # Keep only wanted chromosomes
        if (!is.null(chromList)) {
            message("Removing the following chromosomes: ",
                    paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in% chromList], collapse = ", "))
            SNPset <- SNPset[SNPset$CHROM %in% chromList, ]
        }
        # arrange the chromosomes by natural order sort, eg Chr1, Chr10, Chr2 >>> Chr1, Chr2, Chr10
        SNPset$CHROM <-
            factor(SNPset$CHROM, levels = gtools::mixedsort(unique(SNPset$CHROM)))
        
        # Rename columns
        message("Renaming the following columns: ",
                paste(colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")][grep(highBulk, x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")])], collapse = ", "))
        colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")] <-
            gsub(pattern = highBulk,
                 replacement = "HIGH",
                 x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")])
        
        message("Renaming the following columns: ",
                paste(colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")][grep(lowBulk, x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")])], collapse = ", "))
        colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")] <-
            gsub(pattern = lowBulk,
                 replacement = "LOW",
                 x = colnames(SNPset)[!colnames(SNPset) %in% c("CHROM", "POS", "REF", "ALT")])
        
        # calculate DPs
        SNPset <- SNPset %>%
            dplyr::mutate(
                DP.HIGH = AD_REF.HIGH + AD_ALT.HIGH,
                DP.LOW = AD_REF.LOW + AD_ALT.LOW,
                SNPindex.HIGH = AD_ALT.HIGH / DP.HIGH,
                SNPindex.LOW = AD_ALT.LOW / DP.LOW,
                REF_FRQ = (AD_REF.HIGH + AD_REF.LOW) / (DP.HIGH + DP.LOW),
                deltaSNP = SNPindex.HIGH - SNPindex.LOW
            ) %>%
            dplyr::select(
                -dplyr::contains("HIGH"),-dplyr::contains("LOW"),-dplyr::one_of("deltaSNP", "REF_FRQ"),
                dplyr::matches("AD.*.LOW"),
                dplyr::contains("LOW"),
                dplyr::matches("AD.*.HIGH"),
                dplyr::contains("HIGH"),
                dplyr::everything()
            )
        
        return(as.data.frame(SNPset))
    }


## not exported still only works for GATK...
importFromVCF <- function(file,
                          highBulk = character(),
                          lowBulk = character(),
                          chromList = NULL) {
    
    vcf <- vcfR::read.vcfR(file = file)
    message("Keeping SNPs that pass all filters")
    vcf <- vcf[vcf@fix[, "FILTER"] == "PASS"] 
    
    fix <- dplyr::as_tibble(vcf@fix[, c("CHROM", "POS", "REF", "ALT")]) %>% mutate(Key = seq(1:nrow(.)))
    
    # if (!all(
    #     c(
    #         "CHROM", 
    #         "POS", 
    #         paste0(highBulk, ".AD"), 
    #         paste0(lowBulk, ".AD"), 
    #         paste0(highBulk, ".DP"), 
    #         paste0(lowBulk, ".DP")
    #     ) %in% names(SNPset))) {
    #     stop("One of the required fields is missing. Check your VCF file.")
    # }
    
    tidy_gt <- extract_gt_tidy(vcf, format_fields = c("AD", "DP", "GQ"), gt_column_prepend = "", alleles = FALSE)
    
    SNPset <- tidy_gt %>%
        filter(Indiv == LowBulk) %>% select(-Indiv) %>%
        dplyr::left_join(select(filter(tidy_gt, Indiv == HighBulk),-Indiv),
                         by = "Key",
                         suffix = c(".LOW", ".HIGH")) %>%
        tidyr::separate(
            col = "AD.LOW",
            into = c("AD_REF.LOW", "AD_ALT.LOW"),
            sep = ",",
            extra = "merge",
            convert = TRUE
        ) %>%
        tidyr::separate(
            col = "AD.HIGH",
            into = c("AD_REF.HIGH", "AD_ALT.HIGH"),
            sep = ",",
            extra = "merge", 
            convert = TRUE
        ) %>%
        dplyr::full_join(x = fix, by = "Key") %>%
        dplyr::mutate(
            AD_ALT.HIGH = DP.HIGH - AD_REF.HIGH,
            AD_ALT.LOW = DP.LOW - AD_REF.LOW,
            SNPindex.HIGH = AD_ALT.HIGH / DP.HIGH,
            SNPindex.LOW = AD_ALT.LOW / DP.LOW,
            REF_FRQ = (AD_REF.HIGH + AD_REF.LOW) / (DP.HIGH + DP.LOW),
            deltaSNP = SNPindex.HIGH - SNPindex.LOW
        ) %>%
        select(-Key)
    #Keep only wanted chromosomes
    if (!is.null(chromList)) {
        message("Removing the following chromosomes: ", paste(unique(SNPset$CHROM)[!unique(SNPset$CHROM) %in% chromList], collapse = ", "))
        SNPset <- SNPset[SNPset$CHROM %in% chromList, ]
    }
    as.data.frame(SNPset)
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
    refAlleleFreq,
    filterAroundMedianDepth,
    minTotalDepth,
    maxTotalDepth,
    minSampleDepth,
    minGQ,
    verbose = TRUE) {
    
    org_count <- nrow(SNPset)
    count <- nrow(SNPset)
    
    # Filter by total reference allele frequency
    if (!missing(refAlleleFreq)) {
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
    if (!missing(minGQ)) {
        if (all(c("GQ.LOW", "GQ.HIGH") %in% names(SNPset))) {
        if (verbose) {
            message("Filtering by Genotype Quality: GQ >= ", minGQ)
        }
        SNPset <-
            dplyr::filter(SNPset, GQ.LOW >= minGQ & GQ.HIGH >= minGQ)
        if (verbose) {
            message("...Filtered ", count - nrow(SNPset), " SNPs")
        }
        count <- nrow(SNPset)} 
        else {
            message("GQ columns not found. Skipping...")
        }
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
    return(as.data.frame(SNPset))
}
