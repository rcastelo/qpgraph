# qpgraph: Estimate gene and eQTL networks from high-throughput expression and genotyping assays

[![Bioconductor Time](http://bioconductor.org/shields/years-in-bioc/qpgraph.svg)](http://bioconductor.org/packages/release/bioc/html/qpgraph.html "How long has been qpgraph in a release of Bioconductor")
[![Bioconductor Downloads](http://bioconductor.org/shields/downloads/qpgraph.svg)](http://bioconductor.org/packages/stats/bioc/qpgraph.html "Percentile (top 5/20/50% or 'available') of downloads over the last 6 full months")
[![Bioconductor Commits](http://bioconductor.org/shields/commits/bioc/qpgraph.svg)](http://bioconductor.org/packages/devel/bioc/html/qpgraph.html#svn_source "Average SVN commits (to the devel branch) per month over the last 6 months")
[![Support posts](http://bioconductor.org/shields/posts/qpgraph.svg)](https://support.bioconductor.org/t/qpgraph/ "Bioconductor support site activity on qpgraph, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")

**Current build status**
- `release` [![Bioconductor Availability](http://bioconductor.org/shields/availability/release/qpgraph.svg)](http://bioconductor.org/packages/release/bioc/html/qpgraph.html#archives "Whether qpgraph release is available on all platforms") 
[![Bioconductor Release Build](http://bioconductor.org/shields/build/release/bioc/qpgraph.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/qpgraph/ "Bioconductor release build")
- `development` [![Bioconductor Availability](http://bioconductor.org/shields/availability/devel/qpgraph.svg)](http://bioconductor.org/packages/devel/bioc/html/qpgraph.html#archives "Whether qpgraph devel is available on all platforms") 
[![Bioconductor Devel Build](http://bioconductor.org/shields/build/devel/bioc/qpgraph.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/qpgraph/ "Bioconductor devel build")

## Installation

This is the __development__ version of the R/Bioconductor package qpgraph. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [http://bioconductor.org/packages/qpgraph](http://bioconductor.org/packages/qpgraph) and follow the instructions there to install it.

If you were really looking for this development version, then to install it you
need first to install the version of R that works with Bioc-devel.
See [here](https://www.bioconductor.org/developers/how-to/useDevel/) for more details.

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("qpgraph")
```

Alternatively, you can install it from GitHub using the [devtools](https://github.com/hadley/devtools "devtools") package.

```r
install.packages("devtools")
library(devtools)
install_github("rcastelo/qpgraph")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **qpgraph**
please use the [Bioconductor support site](http://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **qpgraph**
please use the GitHub issues link at the top-right of this page
([https://github.com/rcastelo/qpgraph/issues](https://github.com/rcastelo/qpgraph/issues)).
