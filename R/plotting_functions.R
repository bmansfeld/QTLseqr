#' Plots different paramaters for QTL identification
#'
#' A wrapper for ggplot to plot genome wide distribution of parameters used to
#' identify QTL.
#'
#' @param SNPset a data frame with SNPs and genotype fields as imported by
#'   \code{ImportFromGATK} and after running \code{GetPrimeStats}
#' @param subset a vector of chromosome names for use in quick plotting of
#'   chromosomes of interest. Defaults to
#'   NULL and will plot all chromosomes in the SNPset
#' @param var character. The paramater for plotting. Must be one of: "nSNPs",
#'   "deltaSNP", "Gprime", "negLog10Pval"
#' @param line boolean. If TRUE will plot line graph. If FALSE will plot points.
#'   Plotting points will take more time.
#' @param plotThreshold boolean. Should we plot the False Discovery Rate
#'   threshold (FDR). Only plots line if var is "Gprime" or "negLogPval"
#' @param q numeric. The q-value to use as the FDR threshold. If too low, no
#'   line will be drawn and a warning will be given.
#' @param ... arguments to pass to ggplot2::geom_line or ggplot2::geom_point for
#'   changing colors etc.
#'
#' @return Plots a ggplot graph for all chromosomes or those requested in
#'   \code{subset}. By setting \code{var} to "nSNPs" the distribution of SNPs
#'   used to calculate G' will be plotted. "deltaSNP" will plot a tri-cube
#'   weighted delta SNP-index for each SNP. "Gprime" will plot the tri-cube
#'   weighted G' value. Setting "negLogPval" will plot the -log10 of the p-value
#'   at each SNP. In "Gprime" and "negLogPval" plots, a genome wide FDR threshold of
#'   q can be drawn by setting "plotThreshold" to TRUE. The defualt is a red
#'   line. If you would like to plot a different line we suggest setting
#'   "plotThreshold" to FALSE and manually adding a line using
#'   ggplot2::geom_hline.
#'
#' @examples p <- plotQTLstats(df_filt_6Mb, var = "Gprime", plotThreshold = TRUE, q = 0.01, subset = c("Chr3","Chr4"))
#' @export plotQTLStats

plotQTLStats <-
    function(SNPset,
        subset = NULL,
        var = "nSNPs",
        line = TRUE,
        plotThreshold = FALSE,
        q = 0.05,
        ...) {
        #get fdr threshold by ordering snps by pval then getting the last pval
        #with a qval < q
        fdrT <- getFDRThreshold(SNPset$pvalue, alpha = q)
        logFdrT <- -log10(fdrT)
        GprimeT <- SNPset[which(SNPset$pvalue == fdrT), "Gprime"]
        
        if (length(fdrT) == 0) {
            warning("The q threshold is too low. No line will be drawn")
        }
        
        if (!all(subset %in% unique(SNPset$CHROM))) {
            whichnot <-
                paste(subset[base::which(!subset %in% unique(SNPset$CHROM))], collapse = ', ')
            stop(paste0("The following are not true chromosome names: ", whichnot))
        }
        
        if (!var %in% c("nSNPs", "deltaSNP", "Gprime", "negLog10Pval"))
            stop(
                "Please choose one of the following variables to plot: \"nSNPs\", \"deltaSNP\", \"Gprime\", \"negLog10Pval\""
            )
        
        #don't plot threshold lines in deltaSNPprime or number of SNPs as they are not relevant
        if ((plotThreshold == TRUE &
                var == "deltaSNP") |
                (plotThreshold == TRUE & var == "nSNPs")) {
            message("FDR threshold is not plotted in deltaSNP or nSNPs plots")
            plotThreshold <- FALSE
        }
        SNPset <-
            if (is.null(subset)) {
                SNPset
            } else {
                SNPset[SNPset$CHROM %in% subset, ]
            }
        
        p <- ggplot2::ggplot(data = SNPset) +
            ggplot2::facet_grid( ~ CHROM, scales = "free_x") +
            ggplot2::scale_x_continuous(labels = format_genomic(), name = "Genomic Position (Mb)") +
            ggplot2::theme(plot.margin = ggplot2::margin(
                b = 10,
                l = 20,
                r = 20,
                unit = "pt"
            ))
        
        if (var == "Gprime") {
            threshold <- GprimeT
            p <- p + ggplot2::ylab("G' value")
        }
        
        if (var == "negLog10Pval") {
            threshold <- logFdrT
            p <-
                p + ggplot2::ylab(expression("-" * log[10] * '(p-value)'))
        }
        
        if (var == "nSNPs") {
            p <- p + ggplot2::ylab("Number of SNPs in window")
        }
        
        if (var == "deltaSNP") {
            var <- "tricubeDeltaSNP"
            p <-
                p + ggplot2::ylab(expression(Delta * 'SNP-index')) +
                ggplot2::ylim(-0.55, 0.55) +
                ggplot2::geom_hline(yintercept = 0,
                    color = "black",
                    alpha = 0.4)
        }
        
        if (line) {
            p <-
                p + ggplot2::geom_line(ggplot2::aes_string(x = "POS", y = var), ...)
        }
        
        if (!line) {
            p <- p + ggplot2::geom_point(ggplot2::aes_string(x = "POS", y = var), ...)
        }
        
        if (plotThreshold == TRUE)
            p <-
            p + ggplot2::geom_hline(
                ggplot2::aes_string(yintercept = "threshold"),
                color = "red",
                size = 1,
                alpha = 0.4
            )
        p
        
    }


#' Plots Gprime distribution
#'
#' Plots a ggplot histogram of the distribution of Gprime with a log normal
#' distribution overlay
#'
#' @param SNPset a data frame with SNPs and genotype fields as imported by
#'   \code{ImportFromGATK} and after running \code{GetPrimeStats}
#' @param outlierFilter one of either "deltaSNP" or "Hampel". Method for 
#'   filtering outlier (ie QTL) regions for p-value estimation
#'   
#' @return Plots a ggplot density estimate of the G' value distribution. It will then
#' overlay an estimated log normal distribution with the same mean and variance
#' as the null G' distribution. This will allow to verify if after filtering your G'
#' value appear to be close to log normally and thus can be used to estimate
#' p-values using the non-parametric estimation method described in Magwene et al. (2011). Breifly,
#' using the natural log of Gprime a median absolute deviation (MAD) is
#' calculated. The Gprime set is trimmed to exclude outlier regions (i.e. QTL)
#' based on Hampel's rule. An estimation of the mode of the trimmed set is
#' calculated using the \code{\link[modeest]{mlv}} function from the package modeest. Finally, the mean
#' and variance of the set are estimated using the median and mode are
#' estimated and used to plot the log normal distribution.
#'
#' @examples plotGprimedist(df_filt_6Mb, outlierFilter = "deltaSNP")
#'
#' @seealso \code{\link{GetPvals}} for how p-values are calculated.
#' @export plotGprimeDist

plotGprimeDist <- function(SNPset, outlierFilter = c("deltaSNP", "Hampel"))
{
    if(outlierFilter == "deltaSNP") {
        trimGprime <- SNPset$Gprime[abs(SNPset$deltaSNP) < 0.1]
    } else {
    # Non-parametric estimation of the null distribution of G'
    
    lnGprime <- log(SNPset$Gprime)
    
    # calculate left median absolute deviation for the trimmed G' prime set
    MAD <-
        median(abs(lnGprime[lnGprime <= median(lnGprime)] - median(lnGprime)))
    
    # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
    trimGprime <-
        SNPset$Gprime[lnGprime - median(lnGprime) <= 5.2 * median(MAD)]
    }
    medianTrimGprime <- median(trimGprime)
    
    # estimate the mode of the trimmed G' prime set using the half-sample method
    modeTrimGprime <-
        modeest::mlv(x = trimGprime, bw = 0.5, method = "hsm")$M
    
    muE <- log(medianTrimGprime)
    varE <- abs(muE - log(modeTrimGprime))
    
    #plot Gprime distrubtion
    p <- ggplot2::ggplot(SNPset) +
        ggplot2::xlim(0, max(SNPset$Gprime) + 1) +
        ggplot2::xlab("G' value") +
        ggplot2::geom_density(ggplot2::aes(x = Gprime, color = "Data")) +
        ggplot2::stat_function(
            fun = dlnorm,
            size = 1,
            args = c(meanlog = muE, sdlog = sqrt(varE)),
            aes(color = paste0("Null distribution \n G' ~ lnN(", round(muE, 2), ",",round(varE, 2), ")"))
        ) +
        ggplot2::scale_colour_manual("Distribution", values = c("black", "blue"))# +
        #ggplot2::annotate(x = 10, y = 0.325, geom="text",  
        #    label = paste0("G' ~ lnN(", round(muE, 2), ",",round(varE, 2), ")"),
        #    color = "blue")
    return(p)
}
