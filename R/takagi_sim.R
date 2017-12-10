# QTLseq simulation functions


#' Randomly calculates an alternate allele frequency within a bulk
#'
#' @param n non-negative integer. The number of individuals in each bulk
#' @param pop the population structure. Defaults to "F2" and assumes "RIL" 
#' population otherwise.
#'
#' @return an alternate allele frequency withing the bulk. Used for simulating 
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
        intervals = c(0.005, 0.025, 0.05, 0.95, 0.975, 0.995)) {
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
        
        CI <- data.frame(t(CI))
        names(CI) <- intervals
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