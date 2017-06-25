plotQTL <-
    function(SNPset,
        subset = NULL,
        var = "nSNPs",
        line = TRUE,
        plotThreshold = FALSE,
        q = 0.05,
        ...) {
        #get fdr threshold by ordering snps by pval then getting the last pval
        #with a qval < q
        tmp <- SNPset[order(SNPset$pval, decreasing = F),]
        fdrT <- tmp[sum(tmp$qval <= q), var]

        if (length(fdrT) == 0) {
            warning("The q threshold is too low. No line will be drawn")
        }

        if (!all(subset %in% unique(SNPset$CHROM))) {
            whichnot <- paste(subset[which(!subset %in% unique(SNPset$CHROM))], collapse = ', ')
            stop(paste0("The following are not true chromosome names: ", whichnot))
        }

        if (!var %in% c("nSNPs", "deltaSNP", "Gprime", "negLogPval"))
            stop(
                "Please choose one of the following variables to plot: \"nSNPs\", \"deltaSNP\", \"Gprime\", \"negLogPval\""
            )

        #don't plot threshold lines in deltaSNPprime or number of SNPs as they are not relevant
        if ((plotThreshold == TRUE &
                var == "deltaSNP") | (plotThreshold == TRUE & var == "nSNPs")) {
            message("FDR threshold is not plotted in deltaSNP or nSNPs plots")
            plotThreshold <- FALSE
        }
        SNPset <-
            if (is.null(subset)) {
                SNPset
            } else {
                SNPset[SNPset$CHROM == subset,]
            }

        p <- ggplot2::ggplot(data = SNPset) +
            facet_grid(~ CHROM, scales = "free_x") +
            scale_x_continuous(labels = format_genomic(),
                name = "Genomic Position") +
            theme(plot.margin = margin(
                b = 10,
                l = 20,
                r = 20,
                unit = "pt"
            ))

        if (var == "Gprime") {
            p <- p + ylab("G' value")
        }

        if (var == "negLogPval") {
            p <-
                p + ylab(expression("-" * log[10] * '(p-value)'))
        }

        if (var == "nSNPs") {
            p <- p + ylab("Number of SNPs in window")
        }

        if (var == "deltaSNP") {
            var <- "deltaSNPprime"
            p <- p + ylab(expression(Delta * 'SNP-index')) +
                ylim(-0.55, 0.55) +
                geom_hline(yintercept = 0,
                    color = "black",
                    alpha = 0.4)
        }

        if (line) {
            p <-
                p + geom_line(aes_string(x = "POS", y = var), size = 2, ...)
        }

        if (!line) {
            p <- p + geom_point(aes_string(x = "POS", y = var), ...)
        }

        if (plotThreshold == TRUE)
            p <-
            p + geom_hline(
                yintercept = fdrT,
                color = "red",
                size = 2,
                alpha = 0.4
            )
        p

    }

#' Plots a histogram of the distribution of Gprime with a log normal
#' distribution overlay

plotGprimedist <- function(SNPset, ModeEstMethod = "hsm")
{
    # Non-parametric estimation of the null distribution of G'

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
        modeest::mlv(x = trimGprime, bw = 0.5, method = ModeEstMethod)$M

    muE <- log(medianTrimGprime)
    varE <- abs(muE - log(modeTrimGprime))

    #plot Gprime distrubtion
    p <- ggplot2::ggplot(SNPset) +
        xlim(0, max(SNPset$Gprime) + 1) +
        xlab("G' value") +
        geom_histogram(aes(x = Gprime, y = ..density..), binwidth = 0.5)  +
        stat_function(
            fun = dlnorm,
            size = 1,
            color = 'blue',
            args = c(meanlog = muE, sdlog = sqrt(varE))
        )
    return(p)
}
