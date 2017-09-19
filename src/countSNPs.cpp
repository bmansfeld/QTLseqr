#include <Rcpp.h>
using namespace Rcpp;

//' Count number of SNPs within a sliding window
//' 
//' For each SNP returns how many SNPs are bracketing it within the set window size
//' 
//' @param POS A numeric vector of genomic positions for each SNP
//' @param windowSize The required window size
//' @export countSNPs_cpp
// [[Rcpp::export]]
NumericVector countSNPs_cpp(NumericVector POS, double windowSize) {
    unsigned int nout=POS.size(), i, left=0, right=0;
    NumericVector out(nout);
    
    for( i=0; i < nout; i++ ) {
        while ((right < nout) & (POS[right + 1] <= POS[i] + windowSize / 2))
            right++;
        
        while (POS[left] <= POS[i] - windowSize / 2)
            left++;
        
        out[i] = right - left + 1 ;
    }
    return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
