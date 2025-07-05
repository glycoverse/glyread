# glyread 0.3.0

## Major changes

* Add `read_byonic_pglycoquant()` for reading results from Byonic identification and pGlycoQuant quantification.
* Add `read_msfragger()` for reading results from MSFragger-Glyco.
* `read_pglyco3_pglycoquant()` now performs protein inference automatically to resolve protein ambiguity.
  This means, the variable information tibble now contains "protein", "gene", and "protein_site" columns.
  The previous ambiguous "proteins", "genes", and "protein_sites" columns are removed.

# glyread 0.2.1

## Minor changes

* The "proteins" column of `read_pglyco3_pglycoquant()` now contains protein uniprot accessions, e.g. "P08185".
  When multiple proteins exist, separated with semicolon, e.g. "P08185;P04196".
