#' Return SNPs in significant regions
#'
#' The function takes a SNP set after calculation of p- and q-values and returns 
#' a list containing all SNPs with q-values below a set alpha. Each entry in the list
#' is a SNP set data frame in a contiguous region with adjusted pvalues lower 
#' than the set false discovery rate alpha.
#'
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param alpha the required false discovery rate alpha
#'
#' @export getSigRegions

getSigRegions <- function(SNPset, alpha = 0.05)
{
    if ("qvalue" %in% colnames(SNPset))
    {
        SigRegions <- list()
        for (x in levels(as.factor(SNPset$CHROM))) {
            chr <- as.data.frame(subset(SNPset, CHROM == x))
            
            runs <- S4Vectors::Rle(chr$qvalue <= alpha)
            runvals <- S4Vectors::runValue(runs)
            starts <- S4Vectors::start(runs)
            ends <- S4Vectors::end(runs)
            
            for (i in 1:S4Vectors::nrun(runs)) {
                SigRegions[[length(SigRegions) + 1]] <- if (runvals[i]) {
                    chr[starts[i]:ends[i],]
                }
            }
        }
        return(SigRegions)
    } else {
        stop("Please first run GetPrimeStats or ParGetPrimeStats to ",
            "calculate q-values")
    }
    
}



#' Export a summarized table of QTL
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param alpha the required false discovery rate alpha
#' @param export logical. If TRUE will export a csv table.
#' @param fileName either a character string naming a file or a connection open for writing. "" indicates output to the console.
#'
#' @return Returns a summarized table of QTL identified. The table contains the following columns:
#' \itemize{
#' \item{id - the QTL identification number}
#' \item{chromosome - The chromosome on which the region was identified} 
#' \item{start - the start position on that chromosome, i.e. the position of the first SNP that passes the FDR threshold}
#' \item{end - the end position} 
#' \item{length - the length in basepairs from start to end of the region}
#' \item{nSNPs - the number of SNPs in the region}
#' \item{avgSNPs_Mb - the average number of SNPs/Mb within that region}
#' \item{peakDeltaSNP - the deltaSNP-index value at the peak summit}
#' \item{maxGprime - the max G' score in the region}
#' \item{meanGprime - the average G' score of that region}
#' \item{sdGprime - the standard deviation of G' within the region}
#' \item{AUCaT - the Area Under the Curve but above the Threshold line, an indicator of how significant or wide the peak is}
#' \item{meanPval - the average p-value in the region}
#' \item{meanQval - the average adjusted p-value in the region}
#'}
#' @export getQTLTable
#'
#' @importFrom dplyr %>%

getQTLTable <-
    function(SNPset,
        alpha = 0.05,
        export = FALSE,
        fileName = "QTL.csv")
    {
        QTL <- getSigRegions(SNPset = SNPset, alpha = alpha)
        fdrT <- getFDRThreshold(SNPset$pvalue, alpha = alpha)
        GprimeT <- SNPset[which(SNPset$pvalue == fdrT), "Gprime"]
        merged_QTL <- dplyr::bind_rows(QTL, .id = "id")
        table <- as.data.frame(
            merged_QTL %>%
                dplyr::group_by(id) %>%
                dplyr::summarise(
                    chromosome = unique(as.character(CHROM)),
                    start = min(POS),
                    end = max(POS),
                    length = max(POS) - min(POS),
                    nSNPs = length(POS),
                    avgSNPs_Mb = round(length(POS) / (max(POS) - min(POS)) * 1e6),
                    peakDeltaSNP = ifelse(mean(tricubeDeltaSNP) >= 0, 
                        max(tricubeDeltaSNP), min(tricubeDeltaSNP)),
                    maxGprime = max(Gprime),
                    meanGprime = mean(Gprime),
                    sdGprime = sd(Gprime),
                    AUCaT = sum(diff(POS) * (((head(Gprime, -1) + tail(Gprime, -1)) / 2
                    ) - GprimeT)),
                    meanPval = mean(pvalue),
                    meanQval = mean(qvalue)
                )
        )
        
        
        if (export) {
            write.csv(file = fileName,
                x = table,
                row.names = FALSE)
        }
        return(table)
    }
