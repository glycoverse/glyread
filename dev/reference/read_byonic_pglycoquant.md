# Read Byonic-pGlycoQuant result

If you used Byonic for intact glycopeptide identification, and used
pGlycoQuant for quantification, this is the function for you. It reads
in a pGlycoQuant result file and returns a
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
object. Currently only label-free quantification is supported. Rows
whose protein field contains the literal marker `>Reverse` are excluded.

## Usage

``` r
read_byonic_pglycoquant(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db",
  parse_structure = TRUE,
  multisite = "expand"
)
```

## Arguments

- fp:

  Local path or HTTP(S) URL of the result file.

- sample_info:

  Local path or HTTP(S) URL of the sample information CSV, or a sample
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

- parse_structure:

  Reserved for compatibility. Both formats report compositions.

- multisite:

  For spectra files, expand multisite glycopeptides into site-specific
  rows (`"expand"`) or remove them (`"drop"`).

## Value

An
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
object.

## Which file to use?

Use `Quant.site.list` (recommended) to import the site-specific glycan
abundances reported by pGlycoQuant. Each row describes a protein
accession, protein site and glycan composition. Use `Quant.spectra.list`
when you need peptide sequences or peptide-level site assignments. The
reader detects the format from the column headers. Files from Byonic
result folder are not needed. For instructions on how to use Byonic and
pGlycoQuant, please refer to the manual:
[pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).

## Glycosylation sites

For each PSM in `Quant.spectra.list`, the reader matches the
monoisotopic residue mass of each reported glycan to a bracketed peptide
modification within 0.001 Da. Each peptide position can carry one
reported glycan. All amino acids are eligible. Protein positions are
`Position + peptide_site - 1`, where `Position` is the start of the
peptide in the protein.

The reader retains assignments shared by all maximum one-to-one matches.
Unresolved sites or site-specific compositions are recorded as `NA`,
with one warning summarizing the affected PSMs. If several unresolved
assignments give the same output feature, its PSM intensity is counted
once.

## Multisite glycopeptides

In spectra files, some glycopeptides can have more than one
glycosylation site. By default, they are expanded into multiple rows
with the same quantification value but different `protein_site` and
`glycan_composition`. Set `multisite = "drop"` to remove multisite
glycopeptides instead.

## Sample information

In site files, sample names are the column headers after
`Empty_Separator`. In spectra files, they are taken from
`Intensity(sample)` headers. The `sample` column in `sample_info` must
match these names after any `sample_name_converter` is applied.

## Quantification

Site files supply one abundance per site-specific glycan and sample.
Spectra files are aggregated by summing PSM intensities per glycopeptide
and sample. In both formats, zero abundance is recorded as `NA`.

## Variable information

The following columns could be found in the variable information tibble:

- `peptide`: character, peptide sequence (spectra files only)

- `peptide_site`: integer, site of glycosylation on peptide (spectra
  files only)

- `protein`: character, protein accession

- `protein_site`: integer, site of glycosylation on protein

- `gene`: character, gene name (symbol), when available

- `glycan_composition`:
  [`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
  glycan compositions.

## See also

[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html)
