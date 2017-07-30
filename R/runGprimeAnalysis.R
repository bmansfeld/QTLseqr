#' Identify QTL using a smoothed G statistic
#'
#' A wrapper function that performs the full G prime analysis to identify QTL.
#' 1) The G statistic is culculated by \link{\code{getG}} statisic G is defined by the equation: \deqn{G = 2*\sum_{i=1}^{4}
#' n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2 * \sum n_i * ln(obs(n_i)/exp(n_i))} 
#' Where for each SNP, \eqn{n_i} from i = 1 to 4 corresponds to the reference
#' and alternate allele depths for each bulk, as described in the following
#' table: \tabular{rcc}{ Allele \tab High Bulk \tab Low Bulk \cr Reference \tab
#' \eqn{n_1} \tab \eqn{n_2} \cr Alternate \tab \eqn{n_3} \tab \eqn{n_4} \cr}
#' ...and \eqn{obs(n_i)} are the observed allele depths as described in the data
#' frame. Method 1 calculates the G statistic using expected values assuming
#' read depth is equal for all alleles in both bulks: \deqn{exp(n_1) = ((n_1 +
#' n_2)*(n_1 + n_3))/(n_1 + n_2 + n_3 + n_4)} \deqn{exp(n_2) = ((n_2 + n_1)*(n_2
#' + n_4))/(n_1 + n_2 + n_3 + n_4)} etc...
#' 2) \link{\code{tricubeGStat}} uses local regression to predict a 
#' tricube smoothed version of the G statistc for each SNP. This works as a
#' weighted average across neighboring SNPs that accounts for Linkage
#' disequilibrium (LD) while minizing noise attributed to SNP calling errors. G
#' values for neighboring SNPs within the window are weighted by physical
#' distance from the focal SNP.
#' 3) \link{\code{getPvals}} then estimates p-values for the
#' weighted G' statistic based on the non-parametric estimation method described
#' in \href{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255}{Magwene et al. 2013}. Breifly, using the natural log of Gprime a median 
#' absolute deviation (MAD) is calculated. The Gprime set is trimmed to exclude 
#' outlier regions (i.e. QTL) based on Hampel's rule. An alternate method for
#' filtering out QTL is proposed using absolute delta SNP indeces greater than
#' 0.1 to filter out potential QTL.An estimation of the mode of the trimmed set
#' is calculated using the \code{\link[modeest]{mlv}} function from the package
#' modeest. Finally, the mean and variance of the set are estimated using the
#' median and mode and p-values are estimated from a log normal distribution.
#'
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param WinSize the window size (in base pairs) bracketing each SNP for which 
#'   to calculate the statitics. Magwene et. al recommend a window size of ~25 
#'   cM, but also recommend optionally trying several window sizes to test if 
#'   peaks are over- or undersmoothed.
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for 
#'   filtering outlier (ie QTL) regions for p-value estimation 
#'
#' @return The supplied SNP set tibble after G' analysis. Includes five new columns. 
#' G - The G statistic for each SNP
#' @export runGprimeAnalysis
#'
#' @examples test
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