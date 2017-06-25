plotQTL <-
    function(SNPset,
        subset = NULL,
        var = "nSNPs",
        line = TRUE,
        plotThreshold = FALSE,
        q = 0.05,
        ...){
        tmp <- SNPset[order(SNPset$pval, decreasing = F), ]
        fdrT <- tmp[sum(tmp$qval <= q), var]

        if (length(fdrT) == 0) {
            warning("The q threshold is too low. No line will be drawn")
        }

        if (!var %in% c("nSNPs", "deltaSNP", "Gprime", "negLogPval"))
            stop(
                "Please choose one of the following variables to plot: \"nSNPs\", \"deltaSNP\", \"Gprime\", \"negLogPval\""
            )

        #don't plot threshold lines in deltaSNPprime or number of SNPs as they are not relevant
        if (var == "deltaSNP" | var == "nSNPs")
            plotThreshold <- FALSE

        SNPset <-
            if (is.null(subset)) {
                SNPset
            } else {
                SNPset[SNPset$CHROM == subset, ]
            }

        p <- ggplot(data = SNPset) +
            facet_grid( ~ CHROM, scales = "free_x") +
            scale_x_continuous(labels = format_genomic(),
                name = "Genomic Position") +
            theme(plot.margin = margin(b = 10, l = 20, r = 20, unit = "pt"))

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
                geom_hline(yintercept = 0, color = "black", alpha = 0.4)
        }

        if (line) {
            p <-
                p + geom_line(aes_string(x = "POS", y = var), size = 2, ...)
        }

        if (!line) {
            p <- p + geom_point(aes_string(x = "POS", y = var), ...)
        }

        if (plotThreshold)
            p <-
            p + geom_hline(
                yintercept = fdrT,
                color = "red",
                size = 2,
                alpha = 0.4
            )
        p

    }

