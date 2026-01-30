# Read GlycanFinder result

GlycanFinder is a software for intact glycopeptide identification. This
function reads in the result file and returns a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object. Currently only label-free quantification is supported.

## Usage

``` r
read_glycan_finder(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db",
  parse_structure = TRUE
)
```

## Arguments

- fp:

  File path of the GlycanFinder result file.

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

- orgdb:

  Name of the OrgDb package to use for UniProt to gene symbol
  conversion. Default is "org.Hs.eg.db".

- parse_structure:

  Logical. Whether to parse glycan structures. If `TRUE` (default),
  glycan structures are parsed and included in the `var_info` as
  `glycan_structure` column. If `FALSE`, structure parsing is skipped
  and structure-related columns are removed.

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Which file to use?

You should use the `lfq/lfq.protein-glycopeptides.csv` file from the
GlycanFinder result folder. This file contains the quantification
information for glycopeptides.

## Sample information

The sample information file should be a `csv` file with the first column
named `sample`, and the rest of the columns being sample information.
The `sample` column must match the sample names in the Area columns of
the GlycanFinder result file (e.g., "C_3", "H_3"), although the order
can be different.

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

- `glycan_structure`:
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
  glycan structures (if `parse_structure = TRUE`).

## See also

[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
