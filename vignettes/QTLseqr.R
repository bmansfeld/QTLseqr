## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ----gprimeanalysis-src, eval = FALSE------------------------------------
#  df_filt <- runGprimeAnalysis(df_filt,
#      windowSize = 1e6,
#      outlierFilter = "deltaSNP",
#      filterThreshold = 0.1)

## ----gprimeanalysis-msg, message = TRUE, warning = FALSE, collapse = TRUE, echo = FALSE----
df_filt <- runGprimeAnalysis(df_filt,
    windowSize = 1e6,
    outlierFilter = "deltaSNP",
    filterThreshold = 0.1)

