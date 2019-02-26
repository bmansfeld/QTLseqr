# QTLseqr 0.7.5
## Bugs
* reading files would sometimes guess wrong. Set the default to col_character. This will help with very large read depth importing.

# QTLseqr 0.7.4
## Bugs
* Fixed a compatibility issue with new versions of the `modeest` package. Please note that from now on QTLseqr requires `modeest (> 2.3.2)`

# QTLseqr 0.7.3
## Updates
* Added a `...` for all functions that use tricubed smoothing functions. So that users can easily pass higher maxk values to `raw.locfit`. 
* Added _"A note about window sizes"_ to the vignette.

# QTLseqr 0.7.2
## Updates
* Added `depthDifference` paramater to `filterSNPs` function. This helps filtering SNPs with high absolute differences in read depth between the bulks. 
* `getQTLTable` now also reports the genomic position of the maximum of each peak. 
* Updates to the vignette about filtering SNPs.

# QTLseqr 0.7.1
## Bug fixes
* Corrected a bug in checking for negative bulksizes

# QTLseqr 0.7.0
## Updates
* Added `importFromTable` function to allow users to import from a delimited file.
* Allowed different size bulks in `runQTLseqAnalysis`.
* Updated vignettes and documentation files. 
* Some documentation link fixes

# QTLseqr 0.6.5
## Bug fixes
* Corrected a bug in import that happend when high or low bulks were named with things that looked like CHROM, POS, ALT or REF. Now the function ignores those columns when renaming to HIGH.xx or LOW.xx. Also better definition of column types to force CHROM to be char and POS to be int.

# QTLseqr 0.6.4
## Bug fixes
* Corrected a windowSize that was set to 1e6 instead of the parameter function in 'runQTLanalysis'. This was causing problems in calculating window depth.

# QTLseqr 0.6.3
## Bug fixes
* Corrected a call to global env variables in `importFromGATK`.
* Fixed issues with bulk names that had periods in them.
* Bug fix in export functions using `table` as variable name.

# QTLseqr 0.6.2
## Bug fixes
* Added responsive x axis brakes. Axis brake labels were getting squashed on small chromosomes.

# QTLseqr 0.6.1
## Updates
* added `plotSimulatedThresholds` function to help users visuallize their confisence intervals
* some manual corrections and mods

# QTLseqr 0.6.0
## Updates
* Added QTLseq analysis functionality
* `plotQTLStats` can now plot confidence intervals in $\Delta (SNP\text{-}index)$ plots
* Export functions run faster and allow for detection of QTL in either "Gprime" or "QTLseq" methods
* removed Bioconductor dependency

# QTLseqr 0.5.8
## Updates
* changed $\Delta SNP\text{-}index$ to $\Delta (SNP\text{-}index)$

# QTLseqr 0.5.7
## Updates
* `plotQTLStats` now allows for chromosome facet shape scaling using the 'scaleChroms' parameter

# QTLseqr 0.5.6
## Bug fixes
* Added package `locfit` to Imports. 

# QTLseqr 0.5.5
## Updates
* `plotGprimeDist` null distribution label more accurate.

# QTLseqr 0.5.4

## Updates
* `plotGprimeDist` now plots histograms of filtered and raw data. overlaid with the null dist. Is easier to interpret. 

# QTLseqr 0.5.3

## New features

* Added a `NEWS.md` file to track changes to the package.
* in `plotGprimeDist` plotting now includes density plots of all data, data after QTL filtering and the null-distribution assuming mean and variance of the filtered set.

## Bug fixes
* Corrected a bug in `getFDRThreshold` function that was using regular p-values and not adjusted pvalues to define threshold
