# Read Byonic-Byologic result

If you used Byonic for intact glycopeptide identification, and used
Byologic for quantification, this is the function for you. It reads in a
result file and returns a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object. Currently only label-free quantification is supported.

## Usage

``` r
read_byonic_byologic(
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

  File path of the pGlycoQuant result file.

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

  name of the OrgDb package to use for UniProt to gene symbol
  conversion. Default is "org.Hs.eg.db".

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Which file to use?

Open the .blgc file in the result folder with PMI-Byos. In the "Peptide
List" panel (usually on the bottom right), click "Export content of the
table to a CSV file" button. The exported .csv file is the file you
should use.

## Multisite glycopeptides

Multisite glycopeptides are supported but their `protein_site` will be
set to `NA` since the exact site of glycosylation cannot be determined
unambiguously.

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

## Aggregation

pGlyco3 performs quantification on the PSM level. This level of
information is too detailed for most downstream analyses. This function
aggregate PSMs into glycopeptides through summation. For each
glycopeptide (unique combination of "peptide", "peptide_site",
"protein", "protein_site", "gene", "glycan_composition",
"glycan_structure"), we sum up the quantifications of all PSMs that
belong to this glycopeptide.

## See also

[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
