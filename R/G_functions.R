# G functions - All functions the manipulate the G statistic

GetGStat <- function(SNPset) {
    # Calculates the G statistic using expected values that are equal to half the
    # read depth of each SNP in each bulk
    SNPset$GStat <- apply(SNPset[,5:ncol(SNPset)], 1, function(x) {
        obs <-
            c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
        exp <-
            c(0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]), 0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]))
        2 * sum(obs * log(obs / exp), na.rm = T)
    })
    SNPset

}


GetGStat2 <- function(SNPset) {
    # Version 2 calculates the G statistic using expected values assuming read depth
    # is equal for all alleles in both bulks
    SNPset$GStat2 <- apply(SNPset[,5:ncol(SNPset)], 1, function(x) {
        obs <-
            c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
        exp <- rep((((x["AD_REF.LOW"]) + (x["AD_ALT.LOW"])) * ((x["AD_REF.HIGH"]) + (x["AD_ALT.HIGH"]))) / sum(obs), 4)
        2 * sum(obs * log(obs / exp), na.rm = T)
    })
    SNPset

}


ParGetWeightedStats <- function(VarSet, WinSize = 1e4)
    # For each SNP calculates the G' statistic, a weighted average of G across neighboring SNPs.
    # G' is calculated over a window of WinSize and weighted with a tricube kernel
    # where weights are defined by physical distance away from the focal SNP

{
    # Calculate the number of cores
    n_cores <- detectCores() - 1

    # Initiate cluster
    cl <- makeCluster(n_cores)
    message("Using ",n_cores, " cores")

    # Create empty dataframe to rebuild VarSet
    SW <- data.frame()

    # Calculte half the tricube kernel weights for WinSize
    Dvector <- abs(seq(
        from = 0,
        to = 1,
        length = (WinSize + 1) / 2
    ))
    KNum <- (1 - Dvector ^ 3) ^ 3

    # Calculate G' for each SNP within each chromosome
    for (x in levels(as.factor(VarSet$CHROM))) {

        chr <- as.data.frame(subset(VarSet, CHROM == x))

        message("Calculating G' for Chr: ",x,"...")

        # create sliding window step bins around each SNP
        bin <-
            data.frame(
                start = chr$POS - (WinSize / 2),
                end = chr$POS + (WinSize / 2),
                focal = chr$POS
            )

        SWdata <- parallel::parApply(cl=cl, X=bin, MARGIN=1,FUN= function(y) {
            # the distance from the focal SNP is calculated for each SNP in the window. One (1) is added as an index
            dfromFocal <-
                abs(chr[y["start"] < chr$POS &
                        chr$POS <= y["end"], "POS"] - y["focal"]) + 1

            # A Kernel weight is given to each SNP including the focal SNP
            KNumWeight <- KNum[dfromFocal]
            Kweight <- KNumWeight / sum(KNumWeight)

            # The wighted G stat is calculated by multiplying by the Kernel Weight
            weightedStats <-
                chr[y["start"] < chr$POS & chr$POS <= y["end"], c("GStat","deltaSNP")] * Kweight

            # Calculate G' for the focal SNP
            c(Gprime=sum(weightedStats$GStat),deltaSNPprime=sum(weightedStats$deltaSNP),nSNPs=nrow(weightedStats))

        })
        SWdata<-t(as.data.frame(SWdata))
        chr <- cbind(chr,SWdata)
        SW <- rbind(SW, chr)

    }
    stopCluster(cl)
    SW
}

GetPvals <- function(VarSet) {
    # Non-parametric estimation of the null distribution of G'

    lnGprime <- log(VarSet$Gprime)

    # calculate left median absolute deviation for the trimmed G' prime set
    MAD <-
        median(abs(lnGprime[lnGprime <= median(lnGprime)] - median(lnGprime)))

    # Trim the G prime set to exclude outlier regions (i.e. QTL) using Hampel's rule
    trimGprime <-
        VarSet$Gprime[lnGprime - median(lnGprime) <= 5.2 * median(MAD)]

    medianTrimGprime <- median(trimGprime)

    # estimate the mode of the trimmed G' prime set using the half-sample method
    modeTrimGprime <-
        modeest::mlv(x = trimGprime, bw = 0.5, method = "hsm")$M

    muE <- log(medianTrimGprime)
    varE <- abs(muE - log(modeTrimGprime))

    VarSet$pval <-
        1 - plnorm(q = VarSet$Gprime,
            meanlog = muE,
            sdlog = sqrt(varE))

    VarSet$qval <- p.adjust(p = VarSet$pval, method = "BH")


    return(VarSet)

}


plotGprimedist <- function(VarSet)
{
    #plot Gprime distrubtion
    ggplot2::ggplot(VarSet) + geom_histogram(aes(x = Gprime, y = ..density..), binwidth = 1)  +
        stat_function(
            fun = dlnorm,
            size = 1,
            color = 'blue',
            args = c(meanlog = muE, sdlog = sqrt(varE))
        )
}
