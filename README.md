# qpgraph: Estimate gene and eQTL networks from high-throughput expression and genotyping assays

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/qpgraph.svg)](https://bioconductor.org/packages/release/bioc/html/qpgraph.html "How long has been qpgraph in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/qpgraph.svg)](https://bioconductor.org/packages/stats/bioc/qpgraph.html "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/bioc/qpgraph.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/qpgraph "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Support posts](https://bioconductor.org/shields/posts/qpgraph.svg)](https://support.bioconductor.org/t/qpgraph/ "Support site activity on qpgraph, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")

**Current build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/qpgraph.svg)](https://bioconductor.org/packages/release/bioc/html/qpgraph.html#archives "Whether qpgraph release is available on all platforms") 
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/qpgraph.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/qpgraph/ "Bioconductor release build")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/qpgraph.svg)](https://bioconductor.org/packages/release/bioc/html/qpgraph.html#since "Number of recursive dependencies needed to install package")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/qpgraph.svg)](https://bioconductor.org/packages/devel/bioc/html/qpgraph.html#archives "Whether qpgraph devel is available on all platforms") 
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/qpgraph.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/qpgraph/ "Bioconductor devel build")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/qpgraph.svg)](https://bioconductor.org/packages/devel/bioc/html/qpgraph.html#since "Number of recursive dependencies needed to install package")

## Installation

This is the __development__ version of the R/Bioconductor package qpgraph. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [https://bioconductor.org/packages/qpgraph](https://bioconductor.org/packages/qpgraph) and follow the instructions there to install it.

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
please use the [Bioconductor support site](https://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **qpgraph**
please use the GitHub issues link at the top-right of this page
([https://github.com/rcastelo/qpgraph/issues](https://github.com/rcastelo/qpgraph/issues)).
