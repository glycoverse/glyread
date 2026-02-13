
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyread <a href="https://glycoverse.github.io/glyread/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/glyread)](https://CRAN.R-project.org/package=glyread)
[![R-universe
version](https://glycoverse.r-universe.dev/glyread/badges/version)](https://glycoverse.r-universe.dev/glyread)
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

### Install glycoverse

We recommend installing the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse), which includes
this package and other core glycoverse packages.

### Install glyread alone

If you don’t want to install all glycoverse packages, you can only
install glyread.

You can install the latest release of glyread from
[r-universe](https://glycoverse.r-universe.dev/glyread)
(**recommended**):

``` r
# install.packages("pak")
pak::repo_add(glycoverse = "https://glycoverse.r-universe.dev")
pak::pkg_install("glyread")
```

Or from [GitHub](https://github.com/glycoverse/glyread):

``` r
pak::pkg_install("glycoverse/glyread@*release")
```

Or install the development version (NOT recommended):

``` r
pak::pkg_install("glycoverse/glyread")
```

**Note:** Tips and troubleshooting for the meta-package
[glycoverse](https://github.com/glycoverse/glycoverse) are also
applicable here: [Installation of
glycoverse](https://github.com/glycoverse/glycoverse#installation).

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
