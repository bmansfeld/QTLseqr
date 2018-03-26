#' Return SNPs in significant regions
#'
#' The function takes a SNP set after calculation of p- and q-values or Takagi confidence intervals and returns
#' a list containing all SNPs with q-values or deltaSNP below a set alpha or confidence intervals, respectively. Each entry in the list
#' is a SNP set data frame in a contiguous region with adjusted pvalues lower
#' than the set false discovery rate alpha.
#'
#' @param SNPset Data frame SNP set containing previously filtered SNPs.
#' @param method either "Gprime" or "QTLseq". The method for detecting significant regions.
#' @param alpha numeric. The required false discovery rate alpha for use with \code{method = "Gprime"}
#' @param interval numeric. For use eith \code{method = "QTLseq"} The Takagi based confidence interval requested. This will find the column named "CI_\*\*", where \*\* is the requested interval, i.e. 99.
#'
#' @export getSigRegions

getSigRegions <-
    function(SNPset,
             method = "Gprime",
             alpha = 0.05,
             interval = 99)
    {
        conf <- paste0("CI_", interval)
        
        if (!method %in% c("Gprime", "QTLseq")) {
            stop("method must be either \"Gprime\" or \"QTLseq\"")
        }
        
        if ((method == "Gprime") & !("qvalue" %in% colnames(SNPset))) {
            stop("Please first use runGprimeAnalysis to calculate q-values")
        }
        
        if ((method == "QTLseq") & !(any(names(SNPset) %in% conf))) {
            stop(
                "Cant find the requested confidence interval. Please check that the requested interval exsits or first use runQTLseqAnalysis to calculate confidence intervals"
            )
        }
        
        #QTL <- getSigRegions(SNPset = SNPset, method = method, interval = interval, alpha = alpha)
        fdrT <- getFDRThreshold(SNPset$pvalue, alpha = alpha)
        GprimeT <- SNPset[which(SNPset$pvalue == fdrT), "Gprime"]
        #merged_QTL <- dplyr::bind_rows(QTL, .id = "id")
        SNPset <- SNPset %>%
            dplyr::group_by(CHROM)
        
        if (method == "QTLseq") {
            qtltable <-
                SNPset %>% dplyr::mutate(passThresh = abs(tricubeDeltaSNP) > abs(!!as.name(conf))) %>%
                dplyr::group_by(CHROM, run = {
                    run = rle(passThresh)
                    rep(seq_along(run$lengths), run$lengths)
                }) %>%
                dplyr::filter(passThresh == T) %>% dplyr::ungroup() %>%
                dplyr::group_by(CHROM) %>% dplyr::group_by(CHROM, qtl = {
                    qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
                }) %>%
                #dont need run variable anymore
                dplyr::select(-run,-qtl,-passThresh)
        } else {
            qtltable <- SNPset %>% dplyr::mutate(passThresh = qvalue <= alpha) %>%
                dplyr::group_by(CHROM, run = {
                    run = rle(passThresh)
                    rep(seq_along(run$lengths), run$lengths)
                }) %>%
                dplyr::filter(passThresh == T) %>% dplyr::ungroup() %>%
                dplyr::group_by(CHROM) %>% dplyr::group_by(CHROM, qtl = {
                    qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
                }) %>%
                #dont need run variable anymore
                dplyr::select(-run,-qtl,-passThresh)
        }
        
        qtltable <- as.data.frame(qtltable)
        qtlList <-
            split(qtltable, factor(
                paste(qtltable$CHROM, qtltable$qtl, sep = "_"),
                levels = gtools::mixedsort(unique(
                    paste(qtltable$CHROM, qtltable$qtl, sep = "_")
                ))
            ))
        
        return(qtlList)
    }


#' Export a summarized table of QTL
#' @param SNPset Data frame SNP set containing previously filtered SNPs
#' @param method either "Gprime" or "QTLseq". The method for detecting significant regions.
#' @param alpha numeric. The required false discovery rate alpha for use with \code{method = "Gprime"}
#' @param interval numeric. For use eith \code{method = "QTLseq"} The Takagi based confidence interval requested. This will find the column named "CI_\*\*", where \*\* is the requested interval, i.e. 99.
#' @param export logical. If TRUE will export a csv table.
#' @param fileName either a character string naming a file or a connection open for writing. "" indicates output to the console.
#'
#' @return Returns a summarized table of QTL identified. The table contains the following columns:
#' \itemize{
#' \item{id - the QTL identification number}
#' \item{chromosome - The chromosome on which the region was identified}
#' \item{start - the start position on that chromosome, i.e. the position of the first SNP that passes the FDR threshold}
#' \item{end - the end position}
#' \item{length - the length in basepairs from start to end of the region}
#' \item{nSNPs - the number of SNPs in the region}
#' \item{avgSNPs_Mb - the average number of SNPs/Mb within that region}
#' \item{peakDeltaSNP - the deltaSNP-index value at the peak summit}
#' \item{maxGprime - the max G' score in the region}
#' \item{meanGprime - the average G' score of that region}
#' \item{sdGprime - the standard deviation of G' within the region}
#' \item{AUCaT - the Area Under the Curve but above the Threshold line, an indicator of how significant or wide the peak is}
#' \item{meanPval - the average p-value in the region}
#' \item{meanQval - the average adjusted p-value in the region}
#'}
#' @export getQTLTable
#'
#' @importFrom dplyr %>%

getQTLTable <-
    function(SNPset,
             method = "Gprime",
             alpha = 0.05,
             interval = 99,
             export = FALSE,
             fileName = "QTL.csv")
    {
        conf <- paste0("CI_", interval)
        
        if (!method %in% c("Gprime", "QTLseq")) {
            stop("method must be either \"Gprime\" or \"QTLseq\"")
        }
        
        if ((method == "Gprime") &
            !("qvalue" %in% colnames(SNPset))) {
            stop("Please first use runGprimeAnalysis to calculate q-values")
        }
        
        if ((method == "QTLseq") & !(any(names(SNPset) %in% conf))) {
            stop(
                "Cant find the requested confidence interval. Please check that the requested interval exsits or first use runQTLseqAnalysis to calculate confidence intervals"
            )
        }
        
        #QTL <- getSigRegions(SNPset = SNPset, method = method, interval = interval, alpha = alpha)
        fdrT <- getFDRThreshold(SNPset$pvalue, alpha = alpha)
        GprimeT <- SNPset[which(SNPset$pvalue == fdrT), "Gprime"]
        #merged_QTL <- dplyr::bind_rows(QTL, .id = "id")
        SNPset <- SNPset %>%
            dplyr::group_by(CHROM)
        
        if (method == "QTLseq") {
            qtltable <-
                SNPset %>% dplyr::mutate(passThresh = abs(tricubeDeltaSNP) > abs(!!as.name(conf))) %>%
                dplyr::group_by(CHROM, run = {
                    run = rle(passThresh)
                    rep(seq_along(run$lengths), run$lengths)
                }) %>%
                dplyr::filter(passThresh == T) %>% dplyr::ungroup() %>%
                dplyr::group_by(CHROM) %>% dplyr::group_by(CHROM, qtl = {
                    qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
                }) %>%
                #dont need run variable anymore
                dplyr::select(-run) %>%
                dplyr::summarize(
                    start = min(POS),
                    end = max(POS),
                    length = end - start,
                    nSNPs = length(POS),
                    avgSNPs_Mb = round(length(POS) / (max(POS) - min(POS)) * 1e6),
                    peakDeltaSNP = ifelse(
                        mean(tricubeDeltaSNP) >= 0,
                        max(tricubeDeltaSNP),
                        min(tricubeDeltaSNP)
                    ),
                    avgDeltaSNP = mean(tricubeDeltaSNP)
                )
        } else {
            qtltable <- SNPset %>% dplyr::mutate(passThresh = qvalue <= alpha) %>%
                dplyr::group_by(CHROM, run = {
                    run = rle(passThresh)
                    rep(seq_along(run$lengths), run$lengths)
                }) %>%
                dplyr::filter(passThresh == T) %>% dplyr::ungroup() %>%
                dplyr::group_by(CHROM) %>% dplyr::group_by(CHROM, qtl = {
                    qtl = rep(seq_along(rle(run)$lengths), rle(run)$lengths)
                }) %>%
                #dont need run variable anymore
                dplyr::select(-run) %>%
                dplyr::summarize(
                    start = min(POS),
                    end = max(POS),
                    length = end - start,
                    nSNPs = length(POS),
                    avgSNPs_Mb = round(length(POS) / (max(POS) - min(POS)) * 1e6),
                    peakDeltaSNP = ifelse(
                        mean(tricubeDeltaSNP) >= 0,
                        max(tricubeDeltaSNP),
                        min(tricubeDeltaSNP)
                    ),
                    avgDeltaSNP = mean(tricubeDeltaSNP),
                    #Gprime stuff
                    maxGprime = max(Gprime),
                    meanGprime = mean(Gprime),
                    sdGprime = sd(Gprime),
                    AUCaT = sum(diff(POS) * (((head(Gprime, -1) + tail(Gprime, -1)) / 2
                    ) - GprimeT)),
                    meanPval = mean(pvalue),
                    meanQval = mean(qvalue)
                )
        }
        
        qtltable <- as.data.frame(qtltable)
        
        if (export) {
            write.csv(file = fileName,
                      x = qtltable,
                      row.names = FALSE)
        }
        return(qtltable)
    }
