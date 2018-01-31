
<!-- README.md is generated from README.Rmd. Please edit that file -->
QTLseqr v0.6.2
==============

QTLseqr is an R package for QTL mapping using NGS Bulk Segregant Analysis.

QTLseqr is still under development and is offered with out any guarantee.

### **For more detailed instructions please read the vignette [here](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf)**

Installation
============

<!-- You can install and update QTLseqr by using our [drat](http://dirk.eddelbuettel.com/code/drat.html) repository hosted on our github page: -->
<!-- ```{r drat-install, eval = FALSE} -->
<!-- install.packages("QTLseqr", repos = "http://bmansfeld.github.io/drat") -->
<!-- ``` -->
<!-- OR You can install QTLseqr from github with: -->
You can install QTLseqr from github with:

``` r
# install devtools first to download packages from github
install.packages("devtools")

# use devtools to install QTLseqr
devtools::install_github("bmansfeld/QTLseqr")
```

**Note:** Apart from regular package dependencies, there are some Bioconductor tools that we use as well, as such you will be prompted to install support for Bioconductor, if you haven't already. QTLseqr makes use of C++ to make some tasks significantly faster (like counting SNPs). Because of this, in order to install QTLseqr from github you will be required to install some compiling tools (Rtools and Xcode, for Windows and Mac, respectively).

### For updates read the [NEWS.md](https://github.com/bmansfeld/QTLseqr/blob/master/NEWS.md)

**If you use QTLseqr in published research, please cite:**

> Mansfeld B.N. and Grumet R, QTLseqr: An R package for bulk segregant analysis with next-generation sequencing *bioRxiv* 208140. doi: <https://doi.org/10.1101/208140>

We also recommend citing the paper for the corresponding method you work with.

QTL-seq method:

> Takagi, H., Abe, A., Yoshida, K., Kosugi, S., Natsume, S., Mitsuoka, C., Uemura, A., Utsushi, &gt; H., Tamiru, M., Takuno, S., Innan, H., Cano, L. M., Kamoun, S. and Terauchi, R. (2013), QTL-seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations. *Plant J*, 74: 174–183. <doi:10.1111/tpj.12105>

G prime method:

> Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. *PLOS Computational Biology* 7(11): e1002255. <https://doi.org/10.1371/journal.pcbi.1002255>

Abstract
--------

Next Generation Sequencing Bulk Segregant Analysis (NGS-BSA) is efficient in detecting quantitative trait loci (QTL). Despite the popularity of NGS-BSA and the R statistical platform, no R packages are currently available for NGS-BSA. We present QTLseqr, an R package for NGS-BSA that identifies QTL using two statistical approaches: QTL-seq and G’. These approaches use a simulation method and a tricube smoothed G statistic, respectively, to identify and assess statistical significance of QTL. QTLseqr, can import and filter SNP data, calculate SNP distributions, relative allele frequencies, G’ values, and log10(p-values), enabling identification and plotting of QTL.

Example figure
--------------

![Example figure](https://github.com/bmansfeld/QTLseqr/raw/master/all_plots.png "Example figure")

The package is an R implementation of the analysis described in The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing Magwene PM, Willis JH, Kelly JK (2011) The Statistics of Bulk Segregant Analysis Using Next Generation Sequencing. PLOS Computational Biology 7(11): e1002255. doi: [10.1371/journal.pcbi.1002255](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002255)

It also incorperates Δ(SNP-index) type analysis as described by Takagi, Hiroki, et al. "QTL‐seq: rapid mapping of quantitative trait loci in rice by whole genome resequencing of DNA from two bulked populations." The Plant Journal 74.1 (2013): 174-183. [doi:10.1111/tpj.12105](http://onlinelibrary.wiley.com/doi/10.1111/tpj.12105/full)

Example
-------

### **For more detailed instructions please read the vignette [here](https://github.com/bmansfeld/QTLseqr/raw/master/vignettes/QTLseqr.pdf)**

This is a basic example which shows you how to import and analyze NGS-BSA data.

``` r

#load the package
library("QTLseqr")

#Set sample and file names
HighBulk <- "SRR834931"
LowBulk <- "SRR834927"
file <- "SNPs_from_GATK.table"

#Choose which chromosomes will be included in the analysis (i.e. exclude smaller contigs)
Chroms <- paste0(rep("Chr", 12), 1:12)

#Import SNP data from file
df <-
    importFromGATK(
        file = file,
        highBulk = HighBulk,
        lowBulk = LowBulk,
        chromList = Chroms
     )

#Filter SNPs based on some criteria
df_filt <-
    filterSNPs(
        SNPset = df,
        refAlleleFreq = 0.20,
        minTotalDepth = 100,
        maxTotalDepth = 400,
        minSampleDepth = 40,
        minGQ = 99
    )


#Run G' analysis
df_filt <- runGprimeAnalysis(
    SNPset = df_filt,
    windowSize = 1e6,
    outlierFilter = "deltaSNP")

#Run QTLseq analysis
df_filt <- runQTLseqAnalysis(
    SNPset = df_filt,
    windowSize = 1e6,
    popStruc = "F2",
    replications = 10000,
    intervals = c(95, 99)
)

#Plot
plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.01)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

#export summary CSV
getQTLTable(SNPset = df_filt, alpha = 0.01, export = TRUE, fileName = "my_BSA_QTL.csv")
```
