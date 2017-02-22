GetPvals <- function(VarSet, plotGdist=FALSE){
#Non-parametric estimation of the null distribution of G'
      
      
      lnGprime <- log(VarSet$Gprime)
      
      #calculate left median absolute deviation for the trimmed G' prime set
      MAD <- median(abs(lnGprime[lnGprime <= median(lnGprime)] - median(lnGprime)))
      
      #Trim the G prime set to exclude outlier regions i.e. QTL using Hampel's rule
      trimGprime <- VarSet$Gprime[ lnGprime - median(lnGprime) <= 5.2*median(MAD)]
      
      medianTrimGprime <- median(trimGprime)
      
      #estimate the mode of the trimmed G' prime set using the half-sample method
      modeTrimGprime <- mlv(x=trimGprime, bw=0.5, method = "hsm")$M
      
      muE <- log(medianTrimGprime)
      varE <- abs(muE - log(modeTrimGprime))
      
      VarSet$pval <- 1-plnorm(q= VarSet$Gprime, meanlog = muE, sdlog = sqrt(varE))
      
      # VarSet$Z <- (abs(log(VarSet$Gprime)-  mean(VarSet$Gprime))/sd(VarSet$Gprime))
      
      VarSet$qval <- p.adjust(p = VarSet$pval, method = "BH")
      
      
      if (plotGdist == T) {
      #plot Gprime distrubtion
      p <- ggplot(VarSet) + geom_histogram(aes(x=Gprime,y = ..density..),binwidth = 1)  +
            stat_function(fun = dlnorm, size=1, color='blue',args = c(meanlog = muE, sdlog= sqrt(varE)))
      }
      print(p)
      return(VarSet)
      
}
