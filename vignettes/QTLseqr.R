## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ----quickStart, eval=FALSE----------------------------------------------
#  
#  #load the package
#  library("QTLseqr")
#  
#  #Set sample and file names
#  HighBulk <- "SRR834931"
#  LowBulk <- "SRR834927"
#  file <- "SNPs_from_GATK.table"
#  #Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
#  Chroms <- paste0(rep("Chr", 12), 1:12)
#  
#  #Import SNP data from file
#  df <-
#      importFromGATK(
#          filename = file,
#          highBulk = HighBulk,
#          lowBulk = LowBulk,
#          chromList = Chroms
#       )
#  
#  #Filter SNPs based on some criteria
#  df_filt <-
#      filterSNPs(
#          df,
#          refAlleleFreq = 0.20,
#          minTotalDepth = 100,
#          maxTotalDepth = 400,
#          minSampleDepth = 40,
#          minGQ = 99
#      )
#  
#  
#  #Run G' analysis
#  df_filt <- runGprimeAnalysis(
#      df_filt,
#      windowSize = 1e6,
#      outlierFilter = "deltaSNP")
#  
#  #Plot
#  plotQTLStats(df_filt, var = "deltaSNP", plotThreshold = TRUE, q = 0.01)
#  plotQTLStats(df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)

## ------------------------------------------------------------------------
library("QTLseqr")
rawData <- system.file("extdata", "Yang_et_al_2013.table", package = "QTLseqr", mustWork = TRUE)

## ---- eval = FALSE-------------------------------------------------------
#  rawData <- "C:/PATH/TO/MY/DIR/My_BSA_data.table"

## ------------------------------------------------------------------------
HighBulk <- "SRR834931"
LowBulk <- "SRR834927"
Chroms <- paste0(rep("Chr", 12), 1:12)

## ---- cache = TRUE-------------------------------------------------------
#import data
df <-
    importFromGATK(
        filename = rawData,
        highBulk = HighBulk,
        lowBulk = LowBulk,
        chromList = Chroms
     )

## ------------------------------------------------------------------------
head(df)

## ---- warning = FALSE----------------------------------------------------
library("ggplot2")
ggplot(data = df) + 
    geom_histogram(aes(x = DP.HIGH + DP.LOW)) + 
    xlim(0,1000)


## ---- warning=FALSE------------------------------------------------------
ggplot(data = df) +
    geom_histogram(aes(x = REF_FRQ))

## ----filtSNPs-source, eval = FALSE, message = FALSE----------------------
#  df_filt <-
#      filterSNPs(
#          SNPset = df,
#          refAlleleFreq = 0.20,
#          minTotalDepth = 100,
#          maxTotalDepth = 400,
#          minSampleDepth = 40,
#          minGQ = 99,
#          verbose = TRUE
#      )

## ----filtSNPs-msgs, message = TRUE, warning = FALSE, collapse = TRUE, echo = FALSE----
df_filt <-
    filterSNPs(
        SNPset = df,
        refAlleleFreq = 0.20,
        minTotalDepth = 100,
        maxTotalDepth = 400,
        minSampleDepth = 40,
        minGQ = 99,
        verbose = TRUE
    )

## ----gprimeanalysis-src, eval = FALSE------------------------------------
#  df_filt <- runGprimeAnalysis(df_filt,
#      windowSize = 1e6,
#      outlierFilter = "deltaSNP")

## ----gprimeanalysis-msg, message = TRUE, warning = FALSE, collapse = TRUE, echo = FALSE----
df_filt <- runGprimeAnalysis(df_filt,
    windowSize = 1e6,
    outlierFilter = "deltaSNP")

## ------------------------------------------------------------------------
head(df_filt)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

