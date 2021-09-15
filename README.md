# atena: analysis of transposable elements in R and Bioconductor

[![R-CMD-check-bioc](https://github.com/functionalgenomics/atena/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/functionalgenomics/atena/actions?query=workflow%3AR-CMD-check-bioc)

[![Support posts](https://bioconductor.org/shields/posts/atena.svg)](https://support.bioconductor.org/t/atena/ "Support site activity on atena, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts.")

[![Codecov test coverage](https://codecov.io/gh/functionalgenomics/atena/branch/master/graph/badge.svg)](https://codecov.io/gh/functionalgenomics/atena?branch=master)

This `atena` package provides access from R/Bioconductor to the following
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

To install the __release__ version of this package please go to its package release landing page at [https://bioconductor.org/packages/atena](https://bioconductor.org/packages/atena) and follow the instructions there to install it.

This github repository contains the __development__ version of the R/Bioconductor package `atena`. This version is unstable and should be used only to test new features. 
To install the development version you first need to install the [development version of Bioconductor](https://bioconductor.org/developers/how-to/useDevel) and then type the following line from the R shell:

```r
BiocManager::install("atena", version = "devel")
```

Alternatively, you can install it from this GitHub repo using the [remotes](https://cran.r-project.org/package=remotes) package, as follows:

```
library(remotes)

install_github("functionalgenomics/atena")
```

Provided that you have installed first all its Bioconductor dependencies;
see the `DESCRIPTION` file. The vignette contains an example on how to use it.


## Questions, bug reports and issues

For questions and bug reports regarding the __release__ version of **atena**
please use the [Bioconductor support site](https://support.bioconductor.org "Bioconductor support site").
For bug reports and issues regarding this __development__ version of **atena**
please use the GitHub issues [tab](https://github.com/functionalgenomics/atena/issues) at the top-left of this page.

## Contributing

Contributions to the software codebase of atena are welcome as long as contributors abide to the
terms of the [Bioconductor Contributor Code of Conduct](https://bioconductor.org/about/code-of-conduct).
If you want to contribute to the development of GenomicScores please open an
[issue](https://github.com/functionalgenomics/atena/issues) to start discussing your suggestion or, in case of a
bugfix or a straightforward feature, directly a
[pull request](https://github.com/functionalgenomics/atena/pulls).
