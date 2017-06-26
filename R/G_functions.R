# G functions - All functions the manipulate the G statistic

#' Calculates the G statistic - method 1
#'
#' The function is used by \code{\link{ImportFromGATK}} to calculate the G statisic
#'
#' @param SNPset a data frame with SNPs and genotype fields as imported by ImportFromGATK
#' G is defined by the equation:
#' \deqn{G = 2*\sum_{i=1}^{4} n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2 * \sum n_i * ln(obs(n_i)/exp(n_i))}
#' Where for each SNP, \eqn{n_i} from i = 1 to 4 corresponds to the reference and
#' alternate allele depths for each bulk, as described in the following table:
#' \tabular{rcc}{
#' Allele \tab High Bulk \tab Low Bulk \cr
#' Reference \tab \eqn{n_1} \tab \eqn{n_2} \cr
#' Alternate \tab \eqn{n_3} \tab \eqn{n_4} \cr} ...and \eqn{obs(n_i)} are the observed allele depths as described in the data
#' frame. In this method for calculating G, the expected values \eqn{exp(n_i)}
#' are derived by deviding the read depth for the SNP in each bulk by 2. As we
#' expect 50\% of those reads to support the reference allele.
#' @seealso \code{\link{GetGStat2}}

GetGStat <- function(SNPset) {
    GStat <- apply(SNPset[, 5:ncol(SNPset)], 1, function(x) {
        obs <-
            c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
        exp <-
            c(0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]), 0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]))
        2 * sum(obs * log(obs / exp), na.rm = T)
    })
    GStat
}


#' Calculates the G statistic - method 2
#'
#' The function is used by \code{\link{ImportFromGATK}} to calculate the G statisic
#'
#' @param SNPset a data frame with SNPs and genotype fields as imported by ImportFromGATK
#' G is defined by the equation:
#' \deqn{G = 2*\sum_{i=1}^{4} n_{i}*ln\frac{obs(n_i)}{exp(n_i)}}{G = 2 * \sum n_i * ln(obs(n_i)/exp(n_i))}
#' Where for each SNP, \eqn{n_i} from i = 1 to 4 corresponds to the reference and
#' alternate allele depths for each bulk, as described in the following table:
#' \tabular{rcc}{
#' Allele \tab High Bulk \tab Low Bulk \cr
#' Reference \tab \eqn{n_1} \tab \eqn{n_2} \cr
#' Alternate \tab \eqn{n_3} \tab \eqn{n_4} \cr} ...and \eqn{obs(n_i)} are the
#' observed allele depths as described in the data
#' frame. Method 2 calculates the G statistic using expected values assuming read depth
#' is equal for all alleles in both bulks: \eqn{((n_1 + n_3)*(n_2 + n_4))/(n_1 + n_2 + n_3 + n_4)}
#' @seealso \code{\link{GetGStat}}

GetGStat2 <- function(SNPset) {
    GStat2 <-
        apply(SNPset[, 5:ncol(SNPset)], 1, function(x) {
            obs <-
                c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
            exp <-
                rep((((x["AD_REF.LOW"]) + (x["AD_ALT.LOW"])) * ((x["AD_REF.HIGH"]) + (x["AD_ALT.HIGH"]))) / sum(obs), 4)
            2 * sum(obs * log(obs / exp), na.rm = T)
        })
    GStat2

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

GetPrimeStats <- function(SNPset, WinSize = 1e4, ...)
{
    # Create empty dataframe to rebuild SNPset
    SW <- data.frame()

    # Calculte half the tricube kernel weights for WinSize
    Dvector <- abs(seq(
        from = 0,
        to = 1,
        length = (WinSize + 1) / 2
    ))
    KNum <- (1 - Dvector ^ 3) ^ 3

    # Calculate G' for each SNP within each chromosome
    for (x in levels(as.factor(SNPset$CHROM))) {
        chr <- as.data.frame(subset(SNPset, CHROM == x))

        message("Calculating G' for Chr: ", x, "...")

        # create sliding window step bins around each SNP
        bin <-
            data.frame(
                start = chr$POS - (WinSize / 2),
                end = chr$POS + (WinSize / 2),
                focal = chr$POS
            )

        SWdata <- apply(
            X = bin,
            MARGIN = 1,
            FUN = function(y) {
                # the distance from the focal SNP is calculated for each SNP in the window. One (1) is added as an index
                dfromFocal <-
                    abs(chr[y["start"] < chr$POS &
                            chr$POS <= y["end"], "POS"] - y["focal"]) + 1

                # A Kernel weight is given to each SNP including the focal SNP
                KNumWeight <- KNum[dfromFocal]
                Kweight <- KNumWeight / sum(KNumWeight)

                # The wighted G stat is calculated by multiplying by the Kernel Weight
                weightedStats <-
                    chr[y["start"] < chr$POS &
                            chr$POS <= y["end"], c("GStat", "deltaSNP")] * Kweight

                # Calculate G' for the focal SNP
                c(
                    Gprime = sum(weightedStats$GStat),
                    deltaSNPprime = sum(weightedStats$deltaSNP),
                    nSNPs = nrow(weightedStats)
                )

            }
        )
        SWdata <- t(as.data.frame(SWdata))
        chr <- cbind(chr, SWdata)
        SW <- rbind(SW, chr)

    }
    #calculate p- and q-values
    message("#calculating p- and q-values")
    SW <- GetPvals(SW, ...)
    SW
}


#' Calculate tricube weighted statistics for each SNP - Parallelized version
#'
#' Needs the paralle package to work!
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
#' @param n_cores the number of cores to be used for the process. Defults to
#'   NULL which will then use the number of available cores - 1.
#' @param ... Further arguments passed to \code{modeest::mlv} via
#'   \code{GetPvals}
#' @examples df_filt_4mb <- ParGetPrimeStats(df_filt, WinSize = 4e6, n_cores = 3)
#' @seealso \code{\link{GetPvals}} for how p-values are calculated.
#' @export ParGetPrimeStats

ParGetPrimeStats <- function(SNPset,
    WinSize = 1e6,
    n_cores = NULL,
    ...)

{
    if (!requireNamespace("parallel", quietly = T)) {
        stop(
            "The package 'parallel' needs to be installed for this to run.
            Please install it or use the GetWeightedStats function instead.",
            call. = F
        )
    }
    # Calculate the number of cores
    if (is.null(n_cores)) {
        n_cores <- parallel::detectCores() - 1
    }
    # Initiate cluster
    cl <- parallel::makeCluster(n_cores)
    message("Using ", n_cores, " cores")

    # Create empty dataframe to rebuild SNPset
    SW <- data.frame()

    # Calculte half the tricube kernel weights for WinSize
    Dvector <- abs(seq(
        from = 0,
        to = 1,
        length = (WinSize + 1) / 2
    ))
    KNum <- (1 - Dvector ^ 3) ^ 3

    # Calculate G' for each SNP within each chromosome
    for (x in levels(as.factor(SNPset$CHROM))) {
        chr <- as.data.frame(subset(SNPset, CHROM == x))

        message("Calculating G' for Chr: ", x, "...")

        # create sliding window step bins around each SNP
        bin <-
            data.frame(
                start = chr$POS - (WinSize / 2),
                end = chr$POS + (WinSize / 2),
                focal = chr$POS
            )

        SWdata <-
            parallel::parApply(
                cl = cl,
                X = bin,
                MARGIN = 1,
                FUN = function(y) {
                    # the distance from the focal SNP is calculated for each SNP in the window.
                    # One (1) is added because the focal SNP is stored in the first value of the vector
                    dfromFocal <-
                        abs(chr[y["start"] < chr$POS &
                                chr$POS <= y["end"], "POS"] - y["focal"]) + 1

                    # A Kernel weight is given to each SNP including the focal SNP
                    KNumWeight <- KNum[dfromFocal]
                    Kweight <- KNumWeight / sum(KNumWeight)

                    # The wighted G stat is calculated by multiplying by the Kernel Weight
                    weightedStats <-
                        chr[y["start"] < chr$POS &
                                chr$POS <= y["end"], c("GStat", "deltaSNP")] * Kweight

                    # Calculate G' for the focal SNP
                    c(
                        Gprime = sum(weightedStats$GStat),
                        deltaSNPprime = sum(weightedStats$deltaSNP),
                        nSNPs = nrow(weightedStats)
                    )

                }
            )
        SWdata <- t(as.data.frame(SWdata))
        chr <- cbind(chr, SWdata)
        SW <- rbind(SW, chr)

    }
    parallel::stopCluster(cl)

    #Calculate p- and q-values
    message("Calculating p- and q-values")
    SW <- GetPvals(SW, ...)
    SW
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
#'@export GetPvals

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
    if ("qval" %in% colnames(SNPset))
    {
        SigRegions <- list()
        for (x in levels(as.factor(SNPset$CHROM))) {
            chr <- as.data.frame(subset(SNPset, CHROM == x))

            runs <- S4Vectors::Rle(chr$qval <= alpha)
            runvals <- S4Vectors::runValue(runs)
            starts <- S4Vectors::start(runs)
            ends <- S4Vectors::end(runs)

            for (i in 1:nrun(runs)) {
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
