plotSNPdist <- function(VarSet, WinSize=1e6){
      #Plots the number of SNPs per WinSize
      
      SW<-data.frame()
      for (x in levels(as.factor(VarSet$CHROM))) {
            
            chr<-as.data.frame(subset(VarSet, CHROM == x))
            
            #create sliding window step bins
            bin <- seq(1,max(chr[,"POS"]),WinSize) 
            
            count<-sapply(1:(length(bin)),
                          function(i)
                                nrow(chr[chr$POS>=bin[i]-(0.5*WinSize)  &
                                               chr$POS < bin[i]+(0.5*WinSize),]))
            window<-data.frame(POS=bin, nSNPs=count, CHROM=x)
            SW<-rbind(SW, window)
      }
      ggplot(data=SW) + geom_line(aes(x=POS,y=nSNPs)) + facet_grid(~CHROM, scales="free_x")
}