#' Plots different paramaters for QTL identification
#'
#' A wrapper for ggplot to plot genome wide distribution of parameters used to
#' identify QTL.
#'
#' @param SNPset a data frame with SNPs and genotype fields as imported by
#'   \code{ImportFromGATK} and after running \code{runGprimeAnalysis} or \code{runQTLseqAnalysis}
#' @param subset a vector of chromosome names for use in quick plotting of
#'   chromosomes of interest. Defaults to
#'   NULL and will plot all chromosomes in the SNPset
#' @param var character. The paramater for plotting. Must be one of: "nSNPs",
#'   "deltaSNP", "Gprime", "negLog10Pval"
#' @param scaleChroms boolean. if TRUE (default) then chromosome facets will be 
#'   scaled to relative chromosome sizes. If FALSE all facets will be equal
#'   sizes. This is basically a convenience argument for setting both scales and 
#'   shape as "free_x" in ggplot2::facet_grid.
#' @param line boolean. If TRUE will plot line graph. If FALSE will plot points.
#'   Plotting points will take more time.
#' @param plotThreshold boolean. Should we plot the False Discovery Rate
#'   threshold (FDR). Only plots line if var is "Gprime" or "negLogPval". 
#' @param plotIntervals boolean. Whether or not to plot the two-sided Takagi confidence intervals in "deltaSNP" plots.
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
        scaleChroms = TRUE,
        line = TRUE,
        plotThreshold = FALSE,
        plotIntervals = FALSE,
        q = 0.05,
        ...) {
        
        #get fdr threshold by ordering snps by pval then getting the last pval
        #with a qval < q
        
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
        #if you need to plot threshold get the FDR, but check if there are any values that pass fdr
        
        GprimeT <- 0
        logFdrT <- 0
        
        if (plotThreshold == TRUE) {
            fdrT <- getFDRThreshold(SNPset$pvalue, alpha = q)
            
            if (is.na(fdrT)) {
                warning("The q threshold is too low. No threshold line will be drawn")
                plotThreshold <- FALSE
                
            } else {
                logFdrT <- -log10(fdrT)
                GprimeT <- SNPset[which(SNPset$pvalue == fdrT), "Gprime"]
            }
        }
        
        SNPset <-
            if (is.null(subset)) {
                SNPset
            } else {
                SNPset[SNPset$CHROM %in% subset,]
            }
        
        p <- ggplot2::ggplot(data = SNPset) +
            ggplot2::scale_x_continuous(breaks = seq(from = 0,to = max(SNPset$POS), by = 10^(floor(log10(max(SNPset$POS))))), labels = format_genomic(), name = "Genomic Position (Mb)") +
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
                p + ggplot2::ylab(expression(Delta * '(SNP-index)')) +
                ggplot2::ylim(-0.55, 0.55) +
                ggplot2::geom_hline(yintercept = 0,
                    color = "black",
                    alpha = 0.4)
            if (plotIntervals == TRUE) {
                
                ints_df <-
                     dplyr::select(SNPset, CHROM, POS, dplyr::matches("CI_")) %>% tidyr::gather(key = "Interval", value = "value",-CHROM,-POS)
                
                p <- p + ggplot2::geom_line(data = ints_df, ggplot2::aes(x = POS, y = value, color = Interval)) +
                    ggplot2::geom_line(data = ints_df, ggplot2::aes(
                        x = POS,
                        y = -value,
                        color = Interval
                    ))
            }
        }
        
        if (line) {
            p <-
                p + ggplot2::geom_line(ggplot2::aes_string(x = "POS", y = var), ...)
        }
        
        if (!line) {
            p <-
                p + ggplot2::geom_point(ggplot2::aes_string(x = "POS", y = var), ...)
        }
        
        if (plotThreshold == TRUE)
            p <-
            p + ggplot2::geom_hline(
                ggplot2::aes_string(yintercept = "threshold"),
                color = "red",
                size = 1,
                alpha = 0.4
            )
        
        if (scaleChroms == TRUE) {
           p <- p + ggplot2::facet_grid(~ CHROM, scales = "free_x", space = "free_x")
        } else {
           p <- p + ggplot2::facet_grid(~ CHROM, scales = "free_x")    
        }
        
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
#' @param filterThreshold The absolute delta SNP index to use to filter out
#'   putative QTL (default = 0.1)
#' @param binwidth The binwidth for the histogram. Recomended and default = 0.5
#'
#' @return Plots a ggplot histogram of the G' value distribution. The raw data
#'   as well as the filtered G' values (excluding putatitve QTL) are plotted. It
#'   will then overlay an estimated log normal distribution with the same mean
#'   and variance as the null G' distribution. This will allow to verify if
#'   after filtering your G' value appear to be close to log normally and thus
#'   can be used to estimate p-values using the non-parametric estimation method
#'   described in Magwene et al. (2011). Breifly, using the natural log of
#'   Gprime a median absolute deviation (MAD) is calculated. The Gprime set is
#'   trimmed to exclude outlier regions (i.e. QTL) based on Hampel's rule. An
#'   estimation of the mode of the trimmed set is calculated using the
#'   \code{\link[modeest]{mlv}} function from the package modeest. Finally, the
#'   mean and variance of the set are estimated using the median and mode are
#'   estimated and used to plot the log normal distribution.
#'
#' @examples plotGprimedist(df_filt_6Mb, outlierFilter = "deltaSNP")
#'
#' @seealso \code{\link{getPvals}} for how p-values are calculated.
#' @export plotGprimeDist


plotGprimeDist <-
    function(SNPset,
        outlierFilter = c("deltaSNP", "Hampel"),
        filterThreshold = 0.1,
        binwidth = 0.5)
    {
        if (outlierFilter == "deltaSNP") {
            trim_df <- SNPset[abs(SNPset$deltaSNP) < filterThreshold, ]
            trimGprime <- trim_df$Gprime
        } else {
            # Non-parametric estimation of the null distribution of G'
            
            lnGprime <- log(SNPset$Gprime)
            
            # calculate left median absolute deviation for the trimmed G' prime set
            MAD <-
                median(abs(lnGprime[lnGprime <= median(lnGprime)] - median(lnGprime)))
            
            # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
            trim_df <-
                SNPset[lnGprime - median(lnGprime) <= 5.2 * median(MAD),]
            trimGprime <- trim_df$Gprime
        }
        medianTrimGprime <- median(trimGprime)
        
        # estimate the mode of the trimmed G' prime set using the half-sample method
        modeTrimGprime <-
            modeest::mlv(x = trimGprime, bw = 0.5, method = "hsm")$M
        
        muE <- log(medianTrimGprime)
        varE <- abs(muE - log(modeTrimGprime))
        
        n <- length(trim_df$Gprime)
        bw <- binwidth
        
        #plot Gprime distrubtion
        p <- ggplot2::ggplot(SNPset) +
            ggplot2::xlim(0, max(SNPset$Gprime) + 1) +
            ggplot2::xlab("G' value") +
            ggplot2::geom_histogram(ggplot2::aes(x = Gprime, fill = "Raw Data"), binwidth = bw) +
            ggplot2::geom_histogram(data = trim_df,
                ggplot2::aes(x = Gprime, fill = "After filtering"),
                binwidth = bw) +
            ggplot2::stat_function(
                ggplot2::aes(color = "black"),
                size = 1,
                fun = function(x, mean, sd, n, bw) {
                    dlnorm(x = x,
                        mean = muE,
                        sd = sqrt(varE)) * n * bw
                },
                args = c(
                    mean = muE,
                    sd = sqrt(varE),
                    n = n,
                    bw = bw
                )
            ) +
            
            # ggplot2::stat_function(
            #     fun = dlnorm * n,
            #     size = 1,
            #     args = c(meanlog = muE, sdlog = sqrt(varE)),
            # ggplot2::aes(
            #     color = paste0(
            #         "Null distribution \n G' ~ lnN(",
            #         round(muE, 2),
            #         ",",
            #         round(varE, 2),
        #         ")"
        #     )
        #     )
        # ) +
        ggplot2::scale_fill_discrete(name = "Distribution") +
            ggplot2::scale_colour_manual(name = "Null distribution" , values = "black", labels = as.expression(bquote(~theta["G'"]~" ~ lnN("*.(round(muE, 2))*","*.(round(varE, 2))*")")))  +
            ggplot2::guides(fill = ggplot2::guide_legend(order = 1, reverse = TRUE))
        
        #ggplot2::annotate(x = 10, y = 0.325, geom="text",
        #    label = paste0("G' ~ lnN(", round(muE, 2), ",",round(varE, 2), ")"),
        #    color = "blue")
        return(p)
    }


#' Plots simulation data for QTLseq analysis
#'
# The method for simulating delta SNP-index confidence interval thresholds
#' as described in Takagi et al., (2013). Genotypes are randomly assigned for
#' each indvidual in the bulk, based on the population structure. The total
#' alternative allele frequency in each bulk is calculated at each depth used to simulate
#' delta SNP-indeces, with a user defined number of bootstrapped replication.
#' The requested confidence intervals are then calculated from the bootstraps.
#' This function plots the simulated confidence intervals by the read depth.
#'
#' @param SNPset optional. Either supply your data set to extract read depths from or supply depth vector.
#' @param popStruc the population structure. Defaults to "F2" and assumes "RIL" otherwise.
#' @param bulkSize non-negative integer. The number of individuals in each bulk
#' @param depth optional integer vector. A read depth for which to replicate SNP-index calls. If read depth is defined SNPset will be ignored.
#' @param replications integer. The number of bootstrap replications.
#' @param filter numeric. An optional minimum SNP-index filter
#' @param intervals numeric vector. Confidence intervals supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95\% confidence interval, 2.5\% on each side.
#'
#' @return Plots a deltaSNP by depth plot. Helps if the user wants to know the the delta SNP index needed to pass a certain CI at a specified depth.
#'
#' @export plotSimulatedThresholds
#'
#' @examples plotSimulatedThresholds <- function(SNPset = NULL, popStruc = "F2", bulkSize = 25,   depth = 1:150, replications = 10000, filter = 0.3, intervals = c(95, 99))

plotSimulatedThresholds <-
    function(SNPset = NULL,
             popStruc = "F2",
             bulkSize,
             depth = NULL,
             replications = 10000,
             filter = 0.3,
             intervals = c(95, 99)) {
        
        if (is.null(depth)) {
            if (!is.null(SNPset)) {
                message(
                    "Variable 'depth' not defined, using min and max depth from data: ",
                    min(SNPset$minDP),
                    "-",
                    max(SNPset$minDP)
                )
                depth <- min(SNPset$minDP):max(SNPset$minDP)
            } else {
                stop("No SNPset or depth supplied")
            }
        }
        
        #convert intervals to quantiles
        if (all(intervals >= 1)) {
            message(
                "Returning the following two sided confidence intervals: ",
                paste(intervals, collapse = ", ")
            )
            quantiles <- (100 - intervals) / 200
        } else {
            stop(
                "Convidence intervals ('intervals' paramater) should be supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95% confidence interval, 2.5% on each side."
            )
        }
        
        CI <-
            simulateConfInt(
                popStruc = popStruc,
                bulkSize = bulkSize,
                depth = depth,
                replications = replications,
                filter = filter,
                intervals = quantiles
            )
        
        CI <-
            tidyr::gather(CI, key = "Interval", value = "deltaSNP",-depth)
        
        ggplot2::ggplot(data = CI) + 
            ggplot2::geom_line(ggplot2::aes(x = depth, y = deltaSNP, color = Interval)) +
            ggplot2::geom_line(ggplot2::aes(x = depth, y = -deltaSNP,color = Interval))
    }    