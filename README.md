# qpgraph: Estimate gene and eQTL networks from high-throughput expression and genotyping assays

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/qpgraph.svg)](https://bioconductor.org/packages/release/bioc/html/qpgraph.html "How long has been qpgraph in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/release/qpgraph.svg)](https://bioconductor.org/packages/stats/bioc/qpgraph/ "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months.")
[![Support posts](https://bioconductor.org/shields/posts/qpgraph.svg)](https://support.bioconductor.org/t/qpgraph/ "Support site activity on qpgraph, last 6 months: answered posts/total posts.")
[![R-CMD-check-bioc](https://github.com/rcastelo/qpgraph/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/rcastelo/qpgraph/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/rcastelo/qpgraph/coverage.svg?branch=master)](https://codecov.io/github/rcastelo/qpgraph?branch=master)


**Current build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/qpgraph.svg)](https://bioconductor.org/packages/release/bioc/html/qpgraph.html#archives "Whether qpgraph release is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/qpgraph.svg)](https://bioconductor.org/packages/release/bioc/html/qpgraph.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/qpgraph.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/qpgraph "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/qpgraph.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/qpgraph/ "Bioconductor release build")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/qpgraph.svg)](https://bioconductor.org/packages/devel/bioc/html/qpgraph.html#archives "Whether qpgraph devel is available on all platforms") 
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/qpgraph.svg)](https://bioconductor.org/packages/devel/bioc/html/qpgraph.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/devel/bioc/qpgraph.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/qpgraph "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/qpgraph.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/qpgraph/ "Bioconductor devel build")

## Installation

This is the __development__ version of the R/Bioconductor package qpgraph. This version is unstable and should be used only to test new features. If you are looking for the __release__ version of this package please go to its package release landing page at [https://bioconductor.org/packages/qpgraph](https://bioconductor.org/packages/qpgraph) and follow the instructions there to install it.

If you were really looking for this development version, then to install it you need first to install the [development version of Bioconductor](https://bioconductor.org/developers/how-to/useDevel) and then type the following line from the R shell:


```r
BiocManager::install("qpgraph", version = "devel")
```

Alternatively, you can install it from GitHub using the [remotes](https://github.com/r-lib/remotes "remotes") package.

```r
install.packages("remotes")
library(remotes)
install_github("rcastelo/qpgraph")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **qpgraph**
please use the [Bioconductor support site](https://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **qpgraph**
please use the GitHub issues [tab](https://github.com/rcastelo/qpgraph/issues) at the top-left of this page.

## Contributing

Contributions to the software codebase of qpgraph are welcome as long as contributors abide to the
terms of the [Bioconductor Contributor Code of Conduct](https://bioconductor.org/about/code-of-conduct).
If you want to contribute to the development of qpgraph please open an
[issue](https://github.com/rcastelo/qpgraph/issues) to start discussing your suggestion or, in case of a
bugfix or a straightforward feature, directly a
[pull request](https://github.com/rcastelo/qpgraph/pulls).
