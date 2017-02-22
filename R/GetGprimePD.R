GetGprime<-function(VarSet, WinSize = 1e4)
#     For each SNP calculates the G' statistic, a weighted average of G across neighboring SNPs.
#     G' is calculated over a window of WinSize and weighted with a tricube kernel 
#     where weights are defined by physical distance away from the focal SNP

      {

#     Create empty dataframe to rebuild VarSet
      SW<-data.frame()
      
#     Calculte half the tricube kernel weights for WinSize
      Dvector <- abs(seq(from =0,to=1,length=(WinSize+1)/2))
      KNum <- (1- Dvector^3)^3
      
#     The denominator should be equal to the sum of (1-D^3)^3 for all positions within 
#     the window, up- and down-stream the focal SNP, thus the cumulative sum is counted 
#     twice minus the focal SNP (so as not to count it twice)
##      KDen <- sum(KNum) + sum(KNum) - KNum[1]
##     K<-KNum/KDen
      
#     Calculate G' for each SNP within each chromosome
      for (x in levels(as.factor(VarSet$CHROM))) {
            
            chr<-as.data.frame(subset(VarSet, CHROM == x))
            
            #create sliding window step bins around each  SNP
            bin <- data.frame(start=chr$POS - (WinSize/2), end=chr$POS + (WinSize/2), focal=chr$POS)
            chr$Gprime<-apply(bin, 1, function(y) {
                  
                  # the distance from the focal SNP is calculated for each SNP in the window. One is added as an index 
                  dfromFocal <- abs(chr[y["start"]<chr$POS & chr$POS <= y["end"],"POS"] - y["focal"]  ) + 1
                  
                  #A Kernel weight is given to each SNP including the focal SNP
                  ##Kweight <- K[dfromFocal]
                  KNumWeight <- KNum[dfromFocal]
                  Kweight <- KNumWeight/sum(KNumWeight)
                  
                  #The wighted G stat is calculated by multiplying by the Kernel Weight
                  weightedG <- chr[y["start"]<chr$POS & chr$POS <= y["end"],"GStat"] * Kweight
                  
                  #Calculate G' for the focal SNP
                  sum(weightedG)
                  
            })
            SW<-rbind(SW,chr)
      }
      SW
}