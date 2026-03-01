# Read the result from StrucGP

StrucGP is a software for intact glycopeptide identification. As StrucGP
doesn't support quantification, this function returns an
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with a binary (0/1) expression matrix indicating whether each
glycopeptide was identified in each sample.

## Usage

``` r
read_strucgp(fp, sample_info = NULL, glycan_type = "N", parse_structure = TRUE)
```

## Arguments

- fp:

  File path of the StrucGP result file.

- sample_info:

  File path of the sample information file (csv), or a sample
  information data.frame/tibble. If `NULL` (default), a simple sample
  information tibble will be created.

- glycan_type:

  Glycan type. One of "N", "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", or
  "O-Glc". Default is "N".

- parse_structure:

  Logical. Whether to parse glycan structures. If `TRUE` (default),
  glycan structures are parsed and included in the `var_info` as
  `glycan_structure` column. If `FALSE`, structure parsing is skipped
  and the structure column is removed.

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object with a binary (0/1) expression matrix, where 1 indicates the
glycopeptide was identified in that sample and 0 indicates it was not.

## Variable information

The following columns could be found in the variable information tibble:

- `peptide`: character, peptide sequence

- `protein`: character, protein accession

- `gene`: character, gene name (symbol)

- `protein_site`: integer, site of glycosylation on protein

- `glycan_composition`:
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
  glycan compositions.

- `glycan_structure`:
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
  glycan structures (if `parse_structure = TRUE`).

## Sample information

The sample information file should be a `csv` file with the first column
named `sample`, and the rest of the columns being sample information.
The `sample` column must match the `file_name` column in the StrucGP
result file, although the order can be different.

You can put any useful information in the sample information file.
Recommended columns are:

- `group`: grouping or conditions, e.g. "control" or "tumor", required
  for most downstream analyses

- `batch`: batch information, required for batch effect correction

## See also

[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
