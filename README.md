
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyread <a href="https://glycoverse.github.io/glyread/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/glyread)](https://CRAN.R-project.org/package=glyread)
[![R-CMD-check](https://github.com/glycoverse/glyread/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glycoverse/glyread/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/glycoverse/glyread/graph/badge.svg)](https://app.codecov.io/gh/glycoverse/glyread)
<!-- badges: end -->

The goal of glyread is to read and process quantification results from
common glycomics and glycoproteomics software, such as Byonic, StrucGP,
or pGlyco, and convert them into a `glyexp::expriment()` object (from
the [glyexp](https://github.com/glycoverse/glyexp) package) for further
analysis.

## Installation

You can install the latest release of glyread from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("glycoverse/glyread@*release")
```

Or install the development version:

``` r
remotes::install_github("glycoverse/glyread")
```

## Documentation

-   Getting started:
    [Here](https://glycoverse.github.io/glyread/articles/glyread.html)
-   Reference:
    [Here](https://glycoverse.github.io/glyread/reference/index.html)

## Role in `glycoverse`

`glyread` is the entry point for the `glycoverse` ecosystem. It provides
a unified interface for reading and processing data from various
glycomics and glycoproteomics software, and converting them into a
`glyexp::experiment()` object for further analysis.

## Example

Used pGlyco3 and pGlycoQuant for glycopeptide identification and
quantification? Try `read_pglyco3_pglycoquant()`!

``` r
library(glyread)

exp <- read_pglyco3_pglycoquant(
  "path/to/pGlyco3-pGlycoQuant/result/Quant.spectra.list",
  sample_info = "path/to/sample_info.csv",
  quant_method = "label-free",
  glycan_type = "N"
)
```
