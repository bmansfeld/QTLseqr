format_genomic <- function(...) {
      # Format a vector of numeric values according
      # to the International System of Units.
      # http://en.wikipedia.org/wiki/SI_prefix
      #
      # Based on code by Ben Tupper
      # https://stat.ethz.ch/pipermail/r-help/2012-January/299804.html
      # Args:
      #   ...: Args passed to format()
      #
      # Returns:
      #   A function to format a vector of strings using
      #   SI prefix notation
      #
      
      function(x) {
            limits <- c(1e0,   1e3, 1e6)
            prefix <- c("","Kb","Mb")
            
            # Vector with array indices according to position in intervals
            i <- findInterval(abs(x), limits)
            
            # Set prefix to " " for very small values < 1e-24
            i <- ifelse(i==0, which(limits == 1e0), i)
            
            paste(format(round(x/limits[i], 1),
                         trim=TRUE, scientific=FALSE, ...),
                  prefix[i])
      }
}