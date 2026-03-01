# Read MSFragger-Glyco result

MSFragger-Glyco is a software for glycopeptide identification and
quantification. This function reads in the result file and returns a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Usage

``` r
read_msfragger(
  dp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL
)
```

## Arguments

- dp:

  The directory path of the MSFragger-Glyco result folder.

- sample_info:

  File path of the sample information file (csv), or a sample
  information data.frame/tibble.

- quant_method:

  Quantification method. Either "label-free" or "TMT".

- glycan_type:

  Glycan type. One of "N", "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", or
  "O-Glc". Default is "N".

- sample_name_converter:

  A function to convert sample names from file paths. The function
  should take a character vector of old sample names and return new
  sample names. Note that sample names in `sample_info` should match the
  new names. If NULL, original names are kept.

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Details

This function uses the "psm.tsv" file in each sample folder. Sample
names are extracted from the file paths. They are the parent of each
"psm.tsv" file. For example, "msfragger_result/H1/psm.tsv" will be named
"H1".

## Variable information

The following columns could be found in the variable information tibble:

- `peptide`: character, peptide sequence

- `peptide_site`: integer, site of glycosylation on peptide

- `protein`: character, protein accession

- `protein_site`: integer, site of glycosylation on protein

- `gene`: character, gene name (symbol)

- `glycan_composition`:
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
  glycan compositions.

## See also

[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
