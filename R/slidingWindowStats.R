slidingWindowStats <- function(VarSet, WinSize = 1e6,StepSize=1e4) {

  SW <- data.frame()
  for (x in levels(as.factor(VarSet$CHROM))) {
    chr <- as.data.frame(subset(VarSet, CHROM == x))
    
    # create sliding window step bins
    bin <- seq(1, max(chr[, "POS"]), StepSize)
    
    count <- sapply(1:(length(bin)),
                    function(i)
                      nrow(chr[chr$POS >= bin[i] - (0.5 * WinSize)  &
                                 chr$POS < bin[i] + (0.5 *
                                                       WinSize), ]))
    
    deltaSNPindex <- sapply(1:(length(bin)),
                     function(i)
                       mean(chr[chr$POS >= bin[i] - (0.5 * WinSize)  &
                                  chr$POS < bin[i] + (0.5 *
                                                        WinSize), "deltaSNP"],na.rm = T) )
    
    window <- data.frame(POS = bin,
                         nSNPs = count,
                         deltaSNPindex = deltaSNPindex,
                         CHROM = x)
    SW <- rbind(SW, window)
  }
  return(SW)
}