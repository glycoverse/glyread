# Read glyco-decipher output

Glyco-Decipher is a software for glycopeptide identification and
quantification. This function reads in the result file and returns a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object. Currently only label-free quantification is supported.

## Usage

``` r
read_glyco_decipher(
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

  File path of the pGlyco3 result file.

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

- orgdb:

  name of the OrgDb package to use for UniProt to gene symbol
  conversion. Default is "org.Hs.eg.db".

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Which file to use?

You should use the "site.csv" file in the result folder. This file
contains "Site" and "Glycan" columns, followed by quantification result
for each sample.

## Protein inference and uncertain sites

Glyco-Decipher reports uncertain sites and proteins. This function
automatically performs protein inference using the parsimony method to
find the leader proteins. This ensures each glycopeptide is uniquely
mapped to a single protein, gene, and glycosite. Besides, uncertain
sites on proteins are assigned `NA`.

## Variable information

The following columns could be found in the variable information tibble:

- `protein`: character, protein accession (after protein inference)

- `protein_site`: integer, site of glycosylation on protein (after
  protein inference)

- `gene`: character, gene name (symbol) (after protein inference)

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
