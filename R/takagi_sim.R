# QTLseq simulation functions


#' Randomly calculates an alternate allele frequency within a bulk
#'
#' @param n non-negative integer. The number of individuals in each bulk
#' @param pop the population structure. Defaults to "F2" and assumes "RIL" 
#' population otherwise.
#'
#' @return an alternate allele frequency within the bulk. Used for simulating 
#' SNP-indeces.
#'
simulateAlleleFreq <- function(n, pop = "F2") {
    if (pop == "F2") {
        mean(sample(
            x = c(0, 0.5, 1),
            size = n,
            prob = c(1, 2, 1),
            replace = TRUE
        ))
    } else {
        mean(sample(
            x = c(0, 1),
            size = n,
            prob = c(1, 1),
            replace = TRUE
        ))
    }
}


#' Simulates a delta SNP-index with replication
#'
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param altFreq1 numeric. The alternate allele frequency for bulk A. 
#' @param altFreq2 numeric. The alternate allele frequency for bulk B. 
#' @param replicates integer. The number of bootstrap replications.
#' @param filter numeric. an optional minimum SNP-index filter
#'
#' @return Returns a vector of length replicates delta SNP-indeces 


simulateSNPindex <-
    function(depth,
        altFreq1,
        altFreq2,
        replicates = 10000,
        filter = NULL) {
        SNPindex_H <- rbinom(replicates, size = depth, altFreq1) / depth
        SNPindex_L <-
            rbinom(replicates, size = depth, altFreq2) / depth
        deltaSNP <- SNPindex_H - SNPindex_L
        if (!is.null(filter)) {
            deltaSNP <- deltaSNP[SNPindex_H >= filter | SNPindex_L >= filter]
        }
        deltaSNP
    }


#' Simulation of delta SPP index confidence intervals 
#' 
#' The method for simulating delta SNP-index confidence interval thresholds
#' as described in Takagi et al., (2013). Genotypes are randomly assigned for
#' each indvidual in the bulk, based on the population structure. The total
#' alternative allele frequency in each bulk is calculated at each depth used to simulate 
#' delta SNP-indeces, with a user defined number of bootstrapped replication.
#' The requested confidence intervals are then calculated from the bootstraps.
#'
#' @param popStruc the population structure. Defaults to "F2" and assumes "RIL" 
#' @param bulkSize non-negative integer. The number of individuals in each bulk
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param replications integer. The number of bootstrap replications.
#' @param filter numeric. An optional minimum SNP-index filter
#' @param intervals numeric vector of probabilities with values in [0,1] 
#' corresponding to the requested confidence intervals 
#'
#' @return A data frame of delta SNP-index thresholds corrisponding to the 
#' requested confidence intervals at the user set depths. 
#' @export simulateConfInt
#'
#' @examples CI <-
#' simulateConfInt(
#'    popStruc = "F2",
#'    bulkSize = 50,
#'    depth = 1:100,
#'    intervals = c(0.05, 0.95, 0.025, 0.975, 0.005, 0.995, 0.0025, 0.9975)
     
simulateConfInt <-
    function(popStruc = "F2",
        bulkSize,
        depth = 1:100,
        replications = 10000,
        filter = 0.3,
        intervals = c(0.05, 0.025)) {
        if (popStruc == "F2") {
            message(
                "Assuming bulks selected from F2 population, with ",
                bulkSize,
                " individuals per bulk."
            )
        } else {
            message(
                "Assuming bulks selected from RIL population, with ",
                bulkSize,
                " individuals per bulk."
            )
        }
        
        #makes a vector of possible alt allele frequencies once. this is then sampled for each replicate
        tmp_freq <-
            replicate(n = replications * 10, simulateAlleleFreq(n = bulkSize, pop = popStruc))
        
        message(
            paste0(
                "Simulating ",
                replications,
                " SNPs with reads at each depth: ",
                min(depth),
                "-",
                max(depth)
            )
        )
        message(paste0(
            "Keeping SNPs with >= ",
            filter,
            " SNP-index in both simulated bulks"
        ))
        
        # tmp allele freqs are sampled to produce 'replicate' numbers of probablities. these 
        # are then used as altFreq probs to simulate SNP index values, per bulk.
        CI <- sapply(
            X = depth,
            FUN = function(x)
            {
                quantile(
                    x = simulateSNPindex(
                        depth = x,
                        altFreq1 = sample(
                            x = tmp_freq,
                            size = replications,
                            replace = TRUE
                        ),
                        altFreq2 = sample(
                            x = tmp_freq,
                            size = replications,
                            replace = TRUE
                        ),
                        replicates = replications,
                        filter = filter
                    ),
                    probs = intervals,
                    names = TRUE
                )
            }
        )
        
        CI <- as.data.frame(CI)
        
        if (length(CI) > 1) {
            CI <- data.frame(t(CI))
        }
        
        names(CI) <- paste0("CI_", 100 - (intervals * 200))
        CI <- cbind(depth, CI)
        
        #to long format for easy plotting
        # tidyr::gather(data = CI,
        #     key = interval,
        #     convert = TRUE,
        #     value = SNPindex,-depth) %>%
        #     dplyr::mutate(Confidence = factor(ifelse(
        #         interval > 0.5,
        #         paste0(round((1 - interval) * 200, digits = 1), "%"),
        #         paste0((interval * 200), "%")
        # )))
        CI
    }


#' Calculates delta SNP confidence intervals for QTLseq analysis
#'
#'The method for simulating delta SNP-index confidence interval thresholds
#' as described in Takagi et al., (2013). Genotypes are randomly assigned for
#' each indvidual in the bulk, based on the population structure. The total
#' alternative allele frequency in each bulk is calculated at each depth used to simulate 
#' delta SNP-indeces, with a user defined number of bootstrapped replication.
#' The requested confidence intervals are then calculated from the bootstraps.
#'
#' @param SNPset The data frame imported by \code{ImportFromGATK} 
#' @param popStruc the population structure. Defaults to "F2" and assumes "RIL" otherwise
#' @param bulkSize non-negative integer. The number of individuals in each bulk
#' @param depth integer. A read depth for which to replicate SNP-index calls.
#' @param replications integer. The number of bootstrap replications.
#' @param filter numeric. An optional minimum SNP-index filter
#' @param intervals numeric vector. Confidence intervals supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95\% confidence interval, 2.5\% on each side. 
#'
#' @return A SNPset data frame with delta SNP-index thresholds corrisponding to the 
#' requested confidence intervals matching the tricube smoothed depth at each SNP.
#' @export runQTLseqAnalysis
#'
#' @examples df_filt <- runQTLseqAnalysis(
#' SNPset = df_filt,
#' windowSize = 1e6,
#' popStruc = "F2",
#' replications = 10000,
#' intervals = c(95, 99)
#' )
#' 
#' 
runQTLseqAnalysis <- function(SNPset, windowSize = 1e6,
    popStruc = "F2",
    bulkSize,
    depth = NULL,
    replications = 10000,
    filter = 0.3,
    intervals = c(95, 99)
    ) {
    
    message("Counting SNPs in each window...")
    SNPset <- SNPset %>%
        dplyr::group_by(CHROM) %>%
        dplyr::mutate(nSNPs = countSNPs_cpp(POS = POS, windowSize = windowSize))
    
    message("Calculating tricube smoothed delta SNP index...")
    SNPset <- SNPset %>%
        dplyr::mutate(tricubeDeltaSNP = tricubeStat(POS = POS, Stat = deltaSNP, windowSize))
    
    #convert intervals to quantiles
    if (all(intervals >= 1)) {
        message("Returning the following two sided confidence intervals: ", paste(intervals, collapse = ", "))
        quantiles <- (100 - intervals) / 200
    } else {
        stop(
            "Convidence intervals ('intervals' paramater) should be supplied as two-sided percentiles. i.e. If intervals = '95' will return the two sided 95% confidence interval, 2.5% on each side."
        )
    }
    
    #calculate min depth per snp between bulks
    SNPset <-
        SNPset %>% 
        dplyr::mutate(minDP = pmin(DP.LOW, DP.HIGH))
    
    SNPset <-
        SNPset %>% 
        dplyr::group_by(CHROM) %>% 
        dplyr::mutate(tricubeDP = floor(tricubeStat(POS, minDP, windowSize = 1e6)))
    
    if (is.null(depth)) {
        message(
            "Variable 'depth' not defined, using min and max depth from data: ",
            min(SNPset$minDP),
            "-",
            max(SNPset$minDP)
        )
        depth <- min(SNPset$minDP):max(SNPset$minDP)
    }
    
    #simualte confidence intervals
    CI <-
        simulateConfInt(
            popStruc = popStruc,
            bulkSize = bulkSize,
            depth = depth,
            replications = replications,
            filter = filter,
            intervals = quantiles
        )
    
    
    #match name of column for easier joining of repeat columns
    names(CI)[1] <- "tricubeDP"
    
    #use join as a quick way to match min depth to matching conf intervals.
    SNPset <-
        dplyr::left_join(x = SNPset,
            y = CI #, commented out becuase of above change. need to remove eventually
            # by = c("tricubeDP" = "depth")
            )
    
    as.data.frame(SNPset)
    
}

