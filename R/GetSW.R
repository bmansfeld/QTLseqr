GetVARSlidingWindow <- function(VarSet,
                                WinSize = 1e6,
                                StepSize = 1e4,
                                minSNPs = 10) {
  
    SW <-
      data.frame(
        POS = numeric(0),
        HighSNPindex = numeric(0),
        LowSNPindex = numeric(0),
        deltaSNPindex = numeric(0),
        CHROM = integer(0),
        DEPTH = integer(0),
        nSNPs = integer(0)
      )
    for (x in paste0(rep("Chr", 7), 1:7)) {
      chr <- as.data.frame(subset(VarSet, CHROM == x))
      
      #create sliding window step bins
      bin <- seq(1, max(chr[, "POS"]), StepSize)
      
      SWHighindex <- sapply(1:(length(bin)),
                            function(i)
                              mean(chr[chr$POS >= bin[i] - (0.5 * WinSize) &
                                         chr$POS < bin[i] + (0.5 * WinSize), "SNPindex.HIGH"], na.rm =
                                     T))
      SWLowindex <- sapply(1:(length(bin)),
                           function(i)
                             mean(chr[chr$POS >= bin[i] - (0.5 * WinSize) &
                                        chr$POS < bin[i] + (0.5 * WinSize), "SNPindex.LOW"], na.rm =
                                    T))
      
      SWDeltaindex <- sapply(1:(length(bin)),
                             function(i)
                               mean(chr[chr$POS >= bin[i] - (0.5 * WinSize) &
                                          chr$POS < bin[i] + (0.5 * WinSize), "deltaSNP"], na.rm =
                                      T))
      #get mean SNPindex per step
      # SWindex2<-sapply(1:(length(bin)),
      #                 function(i)
      #                       mean(chr[chr$POS>=bin[i] &
      #                                      chr$POS < bin[i]+ WinSize, "SNPindex"], na.rm=T))
      #get mean depth of window step
      
      SWdepth <- sapply(1:(length(bin)),
                        function(i)
                          mean(chr[chr$POS >= bin[i] - (0.5 * WinSize) &
                                     chr$POS < bin[i] + (0.5 * WinSize), "minDepth"], na.rm =
                                 T))
      # SWdepth2<-sapply(1:(length(bin)),
      #                 function(i)
      #                       mean(chr[chr$POS>=bin[i] &
      #                                      chr$POS < bin[i]+WinSize, "DP"], na.rm=T))
      
      #get count of SNPs in window step
      count <- sapply(1:(length(bin)),
                      function(i)
                        nrow(chr[chr$POS >= bin[i] - (0.5 * WinSize)  &
                                   chr$POS < bin[i] + (0.5 * WinSize), ]))
      
      chrom <-
        data.frame(
          POS = bin,
          HighSNPindex = SWHighindex,
          LowSNPindex = SWLowindex,
          deltaSNPindex = SWDeltaindex,
          CHROM = x,
          DEPTH = SWdepth,
          nSNPs = count
        )
      SW <- rbind(SW, chrom)
    }
  }