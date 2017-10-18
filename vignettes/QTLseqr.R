## ----setup, echo=FALSE, results="hide"-----------------------------------
knitr::opts_chunk$set(tidy=FALSE, cache=TRUE,
                      dev="png",
                      message=FALSE, error=FALSE, warning=TRUE)

## ----import, cache=TRUE--------------------------------------------------
#import data
df <-
    importFromGATK(
        file = rawData,
        highBulk = HighBulk,
        lowBulk = LowBulk,
        chromList = Chroms
     )

