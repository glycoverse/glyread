# glyread

The goal of glyread is to read and process quantification results from
common glycomics and glycoproteomics software, such as Byonic, StrucGP,
or pGlyco, and convert them into a `glyexp::expriment()` object (from
the [glyexp](https://github.com/glycoverse/glyexp) package) for further
analysis.

## Installation

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

## Documentation

- Getting started:
  [Here](https://glycoverse.github.io/glyread/articles/glyread.html)
- Reference:
  [Here](https://glycoverse.github.io/glyread/reference/index.html)

## Role in `glycoverse`

`glyread` is the entry point for the `glycoverse` ecosystem. It provides
a unified interface for reading and processing data from various
glycomics and glycoproteomics software, and converting them into a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object for further analysis.

## Example

Used pGlyco3 and pGlycoQuant for glycopeptide identification and
quantification? Try
[`read_pglyco3_pglycoquant()`](https://glycoverse.github.io/glyread/reference/read_pglyco3_pglycoquant.md)!

``` r
library(glyread)

exp <- read_pglyco3_pglycoquant(
  "path/to/pGlyco3-pGlycoQuant/result/Quant.spectra.list",
  sample_info = "path/to/sample_info.csv",
  quant_method = "label-free",
  glycan_type = "N"
)
```
