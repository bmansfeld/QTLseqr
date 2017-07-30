#' Identify QTL using a smoothed G statistic
#' 
#' A wrapper for all the functions that perform the full G prime analysis to 
#' identify QTL. The following steps are performed:\cr 1) Genome-wide G 
#' statistics are calculated by \code{\link{getG}} \cr 2) G' - A 
#' tricube-smoothed G statistic is predicted by local regression within each 
#' chromosome using \code{\link{tricubeGStat}} \cr 3) P-values are estimated 
#' based using the non-parametric method described by Magwene et al. 2011 with 
#' the function \code{\link{getPvals}} \cr 4) Negative Log10- and 
#' Benjamini-Hochberg adjusted p-values are calculated using 
#' \code{\link[stats]{p.adjust}}
#' 
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param WinSize the window size (in base pairs) bracketing each SNP for which 
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25 
#'   cM, but also recommend optionally trying several window sizes to test if 
#'   peaks are over- or undersmoothed.
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for 
#'   filtering outlier (ie QTL) regions for p-value estimation
#'   
#' @return The supplied SNP set tibble after G' analysis. Includes five new 
#'   columns: \itemize{\item{G - The G statistic for each SNP} \item{Gprime -
#'   The tricube smoothed G statistic based on the supplied window size} 
#'   \item{pvalue - the pvalue at each SNP calculatd by non-parametric
#'   estimation} \item{negLog10Pval - the -Log10(pvalue) supplied for quick
#'   plotting} \item{qvalue - the Benajamini-Hochberg adjusted p-value}}
#' @export runGprimeAnalysis
#'   
#' @examples df_filt <- runGprimeAnalysis(df_filt,windowSize = 2e6,outlierFilter = "deltaSNP")
#' 
runGprimeAnalysis <- function(SNPset, windowSize = 1e6, outlierFilter = "deltaSNP")
{
    SNPset <- SNPset %>%
        group_by(CHROM) %>%
        mutate(
            G = getG(
                LowRef = AD_REF.LOW,
                HighRef = AD_REF.HIGH,
                LowAlt = AD_ALT.LOW,
                HighAlt = AD_ALT.HIGH
            ),
            Gprime = tricubeGStat(POS, G, windowSize)
        ) %>%
        ungroup() %>%
        mutate(
            pvalue = getPvals(deltaSNP, Gprime, outlierFilter), 
            negLog10Pval = -log10(pvalue),
            qvalue = p.adjust(p = pvalue, method = "BH")
        )
    
    return(SNPset)
}