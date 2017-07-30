#Functions for calculating and manipulating the G statistic

#' Calculates the G statistic - method 1
#' 
#' The function is used by \code{\link{runGprimeAnalysis}} to calculate the G
#' statisic G is defined by the equation: \deqn{G = 2*\sum_{i=1}^{4}
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
#' 
#' @param LowRef A vector of the reference allele depth in the low bulk 
#' @param HighRef A vector of the reference allele depth in the high bulk
#' @param LowAlt A vector of the alternate allele depth in the low bulk
#' @param HighAlt A vector of the alternate allele depth in the high bulk
#' 
#' @return A vector of G statistic values with the same length as 
#'   
#' @seealso \href{https://doi.org/10.1371/journal.pcbi.1002255}{The Statistics
#'   of Bulk Segregant Analysis Using Next Generation Sequencing}

getG <- function(LowRef, HighRef, LowAlt, HighAlt)
{
    exp <- c(
        (LowRef+HighRef)*(LowRef+LowAlt)/(LowRef + HighRef + LowAlt + HighAlt),
        (LowRef+HighRef)*(HighRef+HighAlt)/(LowRef + HighRef + LowAlt + HighAlt),
        (LowRef+LowAlt)*(LowAlt+HighAlt)/(LowRef + HighRef + LowAlt + HighAlt),
        (LowAlt+HighAlt)*(HighRef+HighAlt)/(LowRef + HighRef + LowAlt + HighAlt)
    )
    obs <- c(LowRef, HighRef, LowAlt, HighAlt)
    
    G <- 2 * (rowSums(obs * log(matrix(obs, ncol = 4) / matrix(exp, ncol = 4))))
    return(G)
}

#' Calculate tricube weighted statistics for each SNP
#'
#' For each SNP calculates statistics, weighted average of across neighboring
#' SNPs. To account for Linkage disequilibrium (LD) Stats are calculated over a
#' window of WinSize and weighted with a tricube kernel where weights are
#' defined by physical distance away from the focal SNP
#' @return Returns the supplied SNPset data frame with 5 new columns added:
#' @return \code{Gprime} The weighted G statistic caluculted with a tricube
#'   smoothing kernel
#' @return \code{nSNPs} The number of SNPs in the window used to calculate
#'   Gprime
#' @return \code{deltaSNPprime} The weighted delta SNP statistic calculated with
#'   a tricube smooting kernel
#' @return \code{pval} The p-value calculated by Non-parametric estimation of
#'   the null distribution of G'
#' @return \code{qval} The adjusted q-value after Benjamini-Hochberg adjustment
#'
#' @param SNPset a data frame with SNPs and genotype fields as imported by
#'   \code{ImportFromGATK}
#' @param WinSize the window size (in base pairs) bracketing each SNP for which
#'   to calculate the statitics. Magwene et. al recommend a window size of ~24
#'   cM.
#' @param ... Further arguments passed to \code{modeest::mlv} via
#'   \code{GetPvals}
#' @examples df_filt_4mb <- GetPrimeStats(df_filt, WinSize = 4e6)
#' @seealso \code{\link{GetPvals}} for how p-values are calculated.



tricubeGStat <- function(POS, GStat, windowSize = 2e6)
{
    stats::predict(locfit::locfit(GStat ~ locfit::lp(POS, h = windowSize, deg = 0)), POS)
}


#' Non-parametric estimation of the null distribution of G'
#'
#' The function is used by \code{GetPrimeStats} to estimate p-values for the weighted G' statistic based on the
#' non-parametric estimation method described in Magwene et al. 2013. Breifly,
#' using the natural log of Gprime a median absolute deviation (MAD) is
#' calculated. The Gprime set is trimmed to exclude outlier regions (i.e. QTL)
#' based on Hampel's rule. An estimation of the mode of the trimmed set is
#' calculated using the \code{\link[modeest]{mlv}} function from the package modeest. Finally, the mean
#' and variance of the set are estimated using the median and mode and p-values
#' are estimated from a log normal distribution. Adjusted p-values (q-values)
#' are assigned using the \code{\link[stats]{p.adjust}} function.
#'
#' @param SNPset a data frame with SNPs and genotype fields as imported by
#'   \code{ImportFromGATK} and after running \code{GetPrimeStats}
#' @param ModeEstMethod String. The method for estimation of the mode. Passed on to
#' \code{\link[modeest]{mlv}}. The default is half sample method (hsm). See
#' \code{\link[modeest]{mlv}} for other methods.
#' @param ... Further arguments passed to \code{modeest::mlv}
#' 
#' @export GetPvals


GetPvals <- function(SNPset, ModeEstMethod = "hsm", ...) {

    lnGprime <- log(SNPset$Gprime)

    # calculate left median absolute deviation for the trimmed G' prime set
    MAD <-
        median(abs(lnGprime[lnGprime <= median(lnGprime)] - median(lnGprime)))

    # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
    trimGprime <-
        SNPset$Gprime[lnGprime - median(lnGprime) <= 5.2 * median(MAD)]

    medianTrimGprime <- median(trimGprime)

    # estimate the mode of the trimmed G' prime set using the half-sample method
    modeTrimGprime <-
        modeest::mlv(x = trimGprime, bw = 0.5, method = ModeEstMethod, ...)$M

    muE <- log(medianTrimGprime)
    varE <- abs(muE - log(modeTrimGprime))

    SNPset$pval <-
        1 - plnorm(q = SNPset$Gprime,
            meanlog = muE,
            sdlog = sqrt(varE))

    SNPset$negLogPval <- -log10(SNPset$pval)

    SNPset$qval <- p.adjust(p = SNPset$pval, method = "BH")


    return(SNPset)

}

#' Return SNPs in significant regions
#'
#' The function takes a SNP set after calculation of p- and q-values and returns a list
#' containing all SNPs with q-values below a set alpha. Each entry in the list
#' is a SNP set data frame in a contiguous region with
#'
#'@export GetSigRegions

GetSigRegions <- function(SNPset, alpha = 0.05)
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
