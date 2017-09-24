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

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

