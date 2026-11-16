# Read GlyPep-Quant label-free quantification results

Import a GlyPep-Quant CSV produced from Glyco-Decipher or pGlyco3
identifications. `read_glyco_decipher_glypep_quant()` extracts protein
coordinates from Glyco-Decipher glycopeptide labels such as
`Q86XS8@135;`. `read_pglyco3_glypep_quant()` reads pGlyco3 glycopeptide
accessions such as `Q86XS8` and records their unreported protein sites
as `NA`.

## Usage

``` r
read_glyco_decipher_glypep_quant(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db"
)

read_pglyco3_glypep_quant(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db"
)
```

## Arguments

- fp:

  Local path or HTTP(S) URL of the result file.

- sample_info:

  Local path or HTTP(S) URL of the sample information CSV, or a sample
  information data.frame/tibble.

- quant_method:

  Quantification method. Only `"label-free"` is supported.

- glycan_type:

  Glycan type. One of "N", "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", or
  "O-Glc". Default is "N".

- sample_name_converter:

  A function to convert sample names from file paths. The function
  should take a character vector of old sample names and return new
  sample names. Note that sample names in `sample_info` should match the
  new names. If NULL, original names are kept.

- orgdb:

  name of the OrgDb package to use for UniProt to gene symbol
  conversion. Default is "org.Hs.eg.db".

## Value

A
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
with one row per input record. Row metadata include the inferred leader
`protein`, the candidate accessions in `protein_group`, and
`protein_site` (an integer, or `NA` for absent or ambiguous
coordinates). `reported_glycan` retains the original label, including
any mass shift, alongside `glycan_composition`. Site tables retain
`reported_site`. Glycopeptide tables retain the original
`reported_protein` and `reported_peptide` labels and provide the
unmodified sequence in `peptide`.

## Details

Both readers recognize the summary format from its column headers.

- `site.csv` reports a site-specific glycan per row in `Site` and
  `Glycan`.

- `glycopeptide.csv` reports a glycopeptide per row in `Peptide`,
  `Glycan`, and `Protein`.

All remaining columns contain sample abundances. Empty cells and zeros
become `NA`. Rows are retained in file order, including distinct
modified forms of the same peptide or glycan.

## See also

[Published GlyPep-Quant result files
(PXD057799)](https://www.ebi.ac.uk/pride/archive/projects/PXD057799)

## Examples

``` r
fp <- tempfile(fileext = ".csv")
writeLines(c("Peptide,Glycan,Protein,HeLa_01",
             "AAFHNAVAVVIYNNK,Hex(8)HexNAc(2),Q86XS8@135;,831800"), fp)
se <- read_glyco_decipher_glypep_quant(fp)
#> ℹ Reading data
#> ℹ No sample information file passed in. An empty tibble will be used.
#> ℹ Reading data

#> ℹ Package `org.Hs.eg.db` not installed. Skipping gene symbol conversion.
#> ℹ Reading data

#> ✔ Reading data [277ms]
#> 
writeLines(c("Peptide,Glycan,Protein,HeLa_01",
             "AAFHNAVAVVIYNNK,Hex(8)HexNAc(2),Q86XS8,831800"), fp)
se <- read_pglyco3_glypep_quant(fp)
#> ℹ Reading data
#> ℹ No sample information file passed in. An empty tibble will be used.
#> ℹ Reading data

#> ℹ Package `org.Hs.eg.db` not installed. Skipping gene symbol conversion.
#> ℹ Reading data

#> ✔ Reading data [39ms]
#> 
unlink(fp)
```
