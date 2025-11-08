# Read pGlyco3-pGlycoQuant result

If you used pGlyco3 for intact glycopeptide identification, and used
pGlycoQuant for quantification, this is the function for you. It reads
in a pGlycoQuant result file and returns a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object. Currently only label-free quantification is supported.

## Usage

``` r
read_pglyco3_pglycoquant(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  parse_structure = FALSE
)
```

## Arguments

- fp:

  File path of the pGlycoQuant result file.

- sample_info:

  File path of the sample information file (csv), or a sample
  information data.frame/tibble.

- quant_method:

  Quantification method. Either "label-free" or "TMT".

- glycan_type:

  Glycan type. Either "N" or "O". Default is "N".

- sample_name_converter:

  A function to convert sample names from file paths. The function
  should take a character vector of old sample names and return new
  sample names. Note that sample names in `sample_info` should match the
  new names. If NULL, original names are kept.

- parse_structure:

  Logical. Whether to parse glycan structures. If `TRUE`, glycan
  structures are parsed and included in the `var_info` as
  `glycan_structure` column. If `FALSE` (default), structure parsing is
  skipped and structure-related columns are removed.

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Which file to use?

You should use the "Quant.spectra.list" file in the pGlycoQuant result
folder. Files from pGlyco3 result folder are not needed. For
instructions on how to use pGlyco3 and pGlycoQuant, please refer to the
manual:
[pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).

## Variable information

The following columns could be found in the variable information tibble:

- `peptide`: character, peptide sequence

- `peptide_site`: integer, site of glycosylation on peptide

- `protein`: character, protein accession (after protein inference)

- `protein_site`: integer, site of glycosylation on protein (after
  protein inference)

- `gene`: character, gene name (symbol) (after protein inference)

- `glycan_composition`:
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
  glycan compositions.

- `glycan_structure`:
  [`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html),
  glycan structures (if `parse_structure = TRUE`).

## Sample information

The sample information file should be a `csv` file with the first column
named `sample`, and the rest of the columns being sample information.
The `sample` column must match the `RawName` column in the pGlyco3
result file, although the order can be different.

You can put any useful information in the sample information file.
Recommended columns are:

- `group`: grouping or conditions, e.g. "control" or "tumor", required
  for most downstream analyses

- `batch`: batch information, required for batch effect correction

## Protein inference

pGlyco3 reports protein groups. That is, shared glycopeptides are
reported as a group of proteins separated by ";". This function
automatically performs protein inference using the parsimony method to
find the leader proteins. This ensures each glycopeptide is uniquely
mapped to a single protein, gene, and glycosite.

## Aggregation

pGlyco3 performs quantification on the PSM level. This level of
information is too detailed for most downstream analyses. This function
aggregate PSMs into glycopeptides through summation. For each
glycopeptide (unique combination of "peptide", "peptide_site",
"protein", "protein_site", "gene", "glycan_composition",
"glycan_structure"), we sum up the quantifications of all PSMs that
belong to this glycopeptide.

## Glycan structures

pGlyco3 reports a "plausible structure" for each glycan. You can set
`parse_structure = TRUE` to parse these structures into a
"glycan_structure" column as a
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
vector. However, please take caution with these structures, because
pGlyco3 does not have strict quality control on glycan structure
annotations.

## See also

[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
