
# upbm <img src="man/figures/upbm.png" align="right" alt="" width="160"/>

**This package is still under development.** To give it a try, please
use the [latest stable
release](https://github.com/pkimes/upbm/tree/stable-pre-cleanup).

``` r
devtools::install_github("pkimes/upbm", ref = "stable-pre-cleanup")
```

-----

*upbm* provides functions for loading, organizing, and analyzing protein
binding microarray (PBM) data in *R* using standard
[*Bioconductor*](https://bioconductor.org/) classes. The pacakage was
developed with a particular focus on [universal PBMs
(uPBMs)](https://www.ncbi.nlm.nih.gov/pubmed/16998473) and testing for
differential affinity across proteins.

If you have any suggestions on how we can improve the package, [let us
know](https://github.com/pkimes/upbm/issues)\!

## Installation

*upbm* depends on two closely related packages,
[*upbmAux*](http://github.com/pkimes/upbmAux) and
[*upbmData*](https://github.com/pkimes/upbmData), which contain
auxiliary data, e.g. standard uPBM probe designs, and example uPBM data
used throughout examples and vignettes to illustrate the functions in
the package.

``` r
# Install from GitHub along with suggested packages
BiocManager::install("pkimes/upbm")
BiocManager::install("pkimes/upbmAux")
BiocManager::install("pkimes/upbmData")
```

## Usage

With *upbm*, PBM data is first read from GenePix Results (GPR) files and
organized as a *SummarizedExperiment* object with rows corresponding to
probes, and columns corresponding to individual array scans. Several
functions are included for plotting, normalizing, and processing array
data to k-mer level summaries. Additionally, the package includes
functions for testing and identifying both preferentially bound k-mers
for proteins and differentially bound k-mers across proteins. To
performing statistical inference, replicates of the proteins of interest
must be available.

## Related Work

This package complements the [*Universal PBM Analysis
Suite*](http://the_brain.bwh.harvard.edu/PBMAnalysisSuite/index.html),
written in Perl. The data pre-processing steps implemented in *upbm*
largely follow the steps in the *Analysis Suite* software, but have been
reimplemented in R (and some C++). The *Analysis Suite* can be used to
analyze individual samples separately, and provides k-mer level
summaries as non-parametric enrichment scores (E-scores). This differs
from the parameteric mixed-effects model employed in *upbm*, which
requires replicate samples, but allows for statistical hypothesis
testing. The two software suites can be viewed as alternative approaches
to summarizing protein-DNA binding at the k-mer level using PBM data,
but with major differences in the approaches taken to perform inference.
We hope you find *upbm* useful\!
