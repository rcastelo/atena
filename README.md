# atena: analysis of transposable elements in R and Bioconductor

[![Bioconductor Time](https://bioconductor.org/shields/years-in-bioc/atena.svg)](https://bioconductor.org/packages/release/bioc/html/atean.html "How long has been atena in a release of Bioconductor")
[![Bioconductor Downloads](https://bioconductor.org/shields/downloads/release/atena.svg)](https://bioconductor.org/packages/stats/bioc/atena/ "Ranking by number of downloads. A lower number means the package is downloaded more frequently. Determined within a package type (software, experiment, annotation, workflow) and uses the number of distinct IPs for the last 12 months.")
[![Support posts](https://bioconductor.org/shields/posts/atena.svg)](https://support.bioconductor.org/t/atena/ "Support site activity on atena, last 6 months: answered posts/total posts.")
[![R-CMD-check-bioc](https://github.com/rcastelo/atena/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/rcastelo/atena/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/rcastelo/atena/coverage.svg?branch=devel)](https://app.codecov.io/github/rcastelo/atena?branch=devel)
<img align="right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/atena/atena.png" height="200"/>

**Current Bioconductor build status**
- `release` [![Bioconductor Availability](https://bioconductor.org/shields/availability/release/atena.svg)](https://bioconductor.org/packages/release/bioc/html/atena.html#archives "Whether atena release is available on all platforms")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/release/atena.svg)](https://bioconductor.org/packages/release/bioc/html/atena.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/release/bioc/atena.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/atena "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Release Build](https://bioconductor.org/shields/build/release/bioc/atena.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/atena/ "Bioconductor release build")
- `development` [![Bioconductor Availability](https://bioconductor.org/shields/availability/devel/atena.svg)](https://bioconductor.org/packages/devel/bioc/html/atena.html#archives "Whether atena devel is available on all platforms")
[![Bioconductor Dependencies](https://bioconductor.org/shields/dependencies/devel/atena.svg)](https://bioconductor.org/packages/devel/bioc/html/atena.html#since "Number of recursive dependencies needed to install package")
[![Bioconductor Commits](https://bioconductor.org/shields/lastcommit/devel/bioc/atena.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/atena "Time since last commit, possible values: today, < 1 week, < 1 month, < 3 months, since release, before release")
[![Bioconductor Devel Build](https://bioconductor.org/shields/build/devel/bioc/atena.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/atena/ "Bioconductor devel build")

The `atena` package provides access from R/Bioconductor to the following
methods for the quantification of transposable element (TE) expression:

* Jin Y et al. TEtranscripts: a package for including transposable elements
  in differential expression analysis of RNA-seq datasets.
  Bioinformatics. 2015;31(22):3593-3599. DOI:
  [10.1093/bioinformatics/btv422](https://doi.org/10.1093/bioinformatics/btv422).

* Tokuyama M et al. ERVmap analysis reveals genome-wide transcription of human
  endogenous retroviruses. PNAS. 2018;115(50):12565-12572. DOI:
  [10.1073/pnas.1814589115](https://doi.org/10.1073/pnas.1814589115).

* Bendall et al. Telescope: characterization of the retrotranscriptome by
  accurate estimation of transposable element expression.
  PLOS Comp. Biol. 2019;15(9):e1006453. DOI:
  [10.1371/journal.pcbi.1006453](https://doi.org/10.1371/journal.pcbi.1006453).

## Installation

`atena` has been submitted to Bioconductor.

To install the __release__ version of this package please go to its package
release landing page at
[https://bioconductor.org/packages/atena](https://bioconductor.org/packages/atena)
and follow the instructions there to install it.

This github repository contains the __development__ version of the
R/Bioconductor package `atena`. This version is unstable and should be used
only to test new features. To install the development version you first need to
install the
[development version of Bioconductor](https://bioconductor.org/developers/how-to/useDevel)
and then type the following line from the R shell:

```r
BiocManager::install("atena", version="devel")
```

Alternatively, you can directly install it from this GitHub repo as follows:

```
BiocManager::install("rcastelo/atena")
```

## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **atena**
please use the
[Bioconductor support site](https://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **atena**
please use the GitHub issues [tab](https://github.com/rcastelo/atena/issues) at
the top-left of this page.

## Contributing

Contributions to the software codebase of atena are welcome as long as
contributors abide to the terms of the
[Bioconductor Contributor Code of Conduct](https://bioconductor.org/about/code-of-conduct).
If you want to contribute to the development of GenomicScores please open an
[issue](https://github.com/rcastelo/atena/issues) to start discussing your
suggestion or, in case of a bugfix or a straightforward feature, directly a
[pull request](https://github.com/rcastelo/atena/pulls).

## Funding

This software project is supported in part by the Spanish
[Ministry of Science, Innovation and Universities](https://www.ciencia.gob.es).
