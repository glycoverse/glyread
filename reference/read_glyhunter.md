# Read GlyHunter result

This function reads in a
[GlyHunter](https://github.com/fubin1999/glyhunter) result file and
returns a
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Usage

``` r
read_glyhunter(
  fp,
  sample_info = NULL,
  glycan_type = "N",
  sample_name_converter = NULL
)
```

## Arguments

- fp:

  File path of the GlyHunter result file.

- sample_info:

  File path of the sample information file (csv), or a sample
  information data.frame/tibble.

- glycan_type:

  Glycan type. Either "N" or "O". Default is "N".

- sample_name_converter:

  A function to convert sample names from file paths. The function
  should take a character vector of old sample names and return new
  sample names. Note that sample names in `sample_info` should match the
  new names. If NULL, original names are kept.

## Value

An
[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html)
object.

## Which file to use?

Use the "summary_area.csv" file in the GlyHunter result folder. No edit
is necessary.

## Variable information

The following columns could be found in the variable information tibble:

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

## See also

[`glyexp::experiment()`](https://glycoverse.github.io/glyexp/reference/experiment.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
