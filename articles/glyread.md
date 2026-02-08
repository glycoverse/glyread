# Get Started with glyread

Glycomics and glycoproteomics data find their home in `experiment()`
objects from [glyexp](https://github.com/glycoverse/glyexp)—a tidy,
structured format designed specifically for glycobiology workflows.
Working with glycopeptide identification tools like `pGlyco3` or
`MSFragger-Glyco`? You can seamlessly import your results into
`experiment()` objects with just a few lines of code using `glyread`.

## One function, two files — that’s it

Getting started is straightforward. First, pick the function that
matches your identification software. `glyread` currently plays nicely
with these popular tools:

- Byonic (quantification with Byologic):
  [`read_byonic_byologic()`](https://glycoverse.github.io/glyread/reference/read_byonic_byologic.md)
- Byonic (quantification with pGlycoQuant):
  [`read_byonic_pglycoquant()`](https://glycoverse.github.io/glyread/reference/read_byonic_pglycoquant.md)
- Peaks GlycanFinder:
  [`read_glycan_finder()`](https://glycoverse.github.io/glyread/reference/read_glycan_finder.md)
- Glyco-Decipher:
  [`read_glyco_decipher()`](https://glycoverse.github.io/glyread/reference/read_glyco_decipher.md)
- MSFragger-Glyco:
  [`read_msfragger()`](https://glycoverse.github.io/glyread/reference/read_msfragger.md)
- pGlyco3 (built-in quantification):
  [`read_pglyco3()`](https://glycoverse.github.io/glyread/reference/read_pglyco3.md)
- pGlyco3 (quantification with pGlycoQuant):
  [`read_pglyco3_pglycoquant()`](https://glycoverse.github.io/glyread/reference/read_pglyco3_pglycoquant.md)
- StrucGP (no quantification):
  [`read_strucgp()`](https://glycoverse.github.io/glyread/reference/read_strucgp.md)

Next, gather your two input files.

**File 1: Results file** This is the output from your identification
software. Each `read_*()` function expects a specific file format, so
check the function documentation to ensure you’re selecting the right
one.

**File 2: Sample information (CSV)** A simple two-column table that
tells `glyread` about your experimental design:

- `sample`: Sample names as they appear in your results file (order is
  flexible)
- `group`: Experimental conditions, treatments, or groupings (this is
  the recommended column name, but you have some flexibility here)

Pro tip: Quality control samples should be labeled “QC” in the `group`
column — this helps downstream analysis recognize them appropriately.

## Load your data

With your files ready, importing data is a one-liner. Here’s a practical
example: suppose you used pGlyco3 for identification and pGlycoQuant for
quantification, with results in `pglyco3_result.csv` and sample details
in `samples.csv`:

``` r
exp <- read_pglyco3_pglycoquant("pglyco3_result.csv", sample_info = "samples.csv")
```

That’s it — your data is now ready for analysis in a tidy `experiment()`
object.

## What’s next?

Each `read_*()` function has its own quirks and options — file format
variations, optional parameters, and output customizations. Check the
function-specific documentation with
[`?read_pglyco3`](https://glycoverse.github.io/glyread/reference/read_pglyco3.md)
(or whichever function matches your workflow) to fine-tune your data
import.

Happy glyco-analyzing! 🍬
