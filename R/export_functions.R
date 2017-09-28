#' Return SNPs in significant regions
#'
#' The function takes a SNP set after calculation of p- and q-values and returns a list
#' containing all SNPs with q-values below a set alpha. Each entry in the list
#' is a SNP set data frame in a contiguous region with
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
                    chr[starts[i]:ends[i], ]
                }
            }
        }
        return(SigRegions)
    } else {
        stop("Please first run GetPrimeStats or ParGetPrimeStats to ",
            "calculate q-values")
    }
    
}



#' Export
#'
#' @param SNPset
#' @param alpha
#' @param export
#' @param fileName
#'
#' @return
#' @export getQTLTable
#'
#' @examples
getQTLTable <-
    function(SNPset,
        alpha = 0.05,
        export = FALSE,
        fileName = "QTL.csv")
    {
        QTL <- getSigRegions(SNPset = SNPset, alpha = alpha)
        fdrT <- getFDRThreshold(SNPset$pvalue, alpha = alpha)
        GprimeT <- SNPset[which(SNPset$pvalue == fdrT), "Gprime"]
        table <- as.data.frame(t(
            sapply(
                X = QTL,
                FUN = dplyr::summarise,
                Chromosome = unique(as.character(CHROM)),
                start = min(POS),
                end = max(POS),
                length = max(POS) - min(POS),
                nSNPs = length(POS),
                avgSNPs_Mb = round(length(POS) / (max(POS) - min(POS)) * 1e6),
                meanGprime = mean(Gprime),
                sdGprime = sd(Gprime),
                AUCaT = sum(diff(POS) * (((head(Gprime,-1) + tail(Gprime,-1)) / 2
                ) - GprimeT)),
                meanPval = mean(pvalue),
                meanQval = mean(qvalue)
            )
        ))
        if (export) {
            write.csv(file = fileName,
                x = table,
                row.names = FALSE)
        }
        return(table)
    }
