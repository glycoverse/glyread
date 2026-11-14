# Read pGlyco3-pGlycoQuant result

If you used pGlyco3 for intact glycopeptide identification, and used
pGlycoQuant for quantification, this is the function for you. It reads
in a pGlycoQuant result file and returns a
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
object. Currently only label-free quantification is supported.

## Usage

``` r
read_pglyco3_pglycoquant(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  parse_structure = FALSE,
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

  Quantification method. Either "label-free" or "TMT".

- glycan_type:

  Glycan type. One of "N", "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", or
  "O-Glc". Default is "N".

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

- orgdb:

  Name of the OrgDb package for mapping site-file protein accessions to
  gene symbols. Spectra files supply their own gene annotations.

## Value

An
[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html)
object.

## Which file to use?

Use `Quant.site.list` (recommended) for the reported site-specific
glycan abundances, or `Quant.spectra.list` for peptide-level data. The
reader detects the format from the column headers. Site files report
protein, protein site and glycan composition directly. Protein inference
and PSM aggregation apply only to spectra files.

`Quant.site.list` omits glycan structures. If you need the structures
reported by pGlyco3, read `Quant.spectra.list` with
`parse_structure = TRUE`.

Files from pGlyco3 result folder are not needed. For instructions on how
to use pGlyco3 and pGlycoQuant, please refer to the manual:
[pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).

## Variable information

Site files provide `protein`, `protein_site`, `glycan_composition`, and
`gene` when available. Spectra files can provide the following columns:

- `peptide`: character, peptide sequence

- `peptide_site`: integer, site of glycosylation on peptide

- `protein`: character, protein accession (after protein inference)

- `protein_group`: character, all protein accessions before inference,
  separated by `;`, retained in `rowData`. Accessions are combined
  across records and samples when they map to the same glycopeptide.

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

In site files, sample names are the column headers after
`Empty_Separator`. In spectra files, they are taken from
`Intensity(sample)` headers. The `sample` column in `sample_info` must
match these names after any `sample_name_converter` is applied.

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

[`glyexp::GlycoproteomicSE()`](https://glycoverse.github.io/glyexp/reference/GlycoproteomicSE.html),
[`glyrepr::glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.html),
[`glyrepr::glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.html)
