# glyread (development version)

# glyread 0.8.3

## Minor improvements

* `glycan_type` parameter now supports "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", and "O-Glc". "O" is not a valid value anymore. This is to be in consistent with the update introduced in `glyexp` 0.10.2.

# glyread 0.8.2

## Minor improvements and fixes

* glyread now depends on the CRAN version of glyparse.

# glyread 0.8.1

## Minor improvements and fixes

* glyread now depends on the CRAN version of glyrepr.

# glyread 0.8.0

## New features

* `read_strucgp()` now has a 0/1 expression matrix, indicating whether each glycopeptide was identified in each sample.
* `read_strucgp()` now performs protein inference automatically to resolve protein ambiguity.

## Minor improvements and bug fixes

* Protein inference algorithm is optimized for performance.

# glyread 0.7.0

## Breaking changes

* `read_strucgp()` has been redesigned with a new API. It now returns an [glyexp::experiment()] object instead of a tibble, to be consistent with other `read_xxx()` functions. As StrucGP doesn't support quantification, the expression matrix is filled with `NA` values.

# glyread 0.6.3

## Minor improvements and bug fixes

* Fix the bug in `read_strucgp()` that glycan compositions with "NeuGc" could not be parsed.
* Fix the bug in `read_strucgp()` that glycan structures with additional modifications (e.g. "+Ammonium(+17)") could not be parsed.
* Fix the bug in `read_strucgp()` that sample information was not included in result tibble.

# glyread 0.6.2

## Minor improvements and bug fixes

* Fix bugs introduced by the breaking changes in `glyexp` 0.10.0.

# glyread 0.6.1

## Minor improvements and bug fixes

* Update dependencies to depend on release versions of glycoverse packages.

# glyread 0.6.0

## New features

* Add `read_strucgp()` to read results from StrucGP.
* Add `read_glyco_decipher()` to read results from Glyco-Decipher.

## Minor improvements

* `read_byonic_byologic()` and `read_byonic_pglycoquant()` do not drop glycopeptides with uncertain sites anymore.
  Instead, they fill the "protein_site" column with `NA` for these glycopeptides.

# glyread 0.5.0

## Breaking changes

* Add a `parse_structure` parameter to `read_pglyco3()` and `read_pglyco3_pglycoquant()`. If `FALSE` (default), glycan structure strings are not parsed, and "glycan_structure" column will not be in the variable information tibble. This improves performance when structure-based analysis is not needed, and avoids implicitly implying that the structures given by pGlyco3 is accurate. Well, they are not.

## Minor improvements and bug fixes

* `read_pglyco3()` and `read_pglyco3_pglycoquant()` now correctly parse "pH" and "aH" monosaccharides in glycan compositions and glycan structures.
* Messages about protein inference changes from "Performing protein inference" to "Finding leader proteins" to better illustrate what is actually happening. Relavent documentations are also updated.
* Add a "Varaible information" section to documentations of all functions. This includes the descriptions of columns in the variable information tibble.
* Fix a bug that Uniprot isomers ("O75882-2") could not be parsed correctly.

# glyread 0.4.2

## Minor improvements

* Remove redundant emoji in the "Getting Started" vignette.

# glyread 0.4.1

## Minor improvements

* Add a "Getting Started" vignette.

# glyread 0.4.0

## Major changes

This is a big update! In addition to several new workflows supported,
we redesign the API to make `read_xxx()` functions behave more consistently.
Now the `read_xxx()` functions will try to return experiments with similar variable information format.
To ensure this, some functions perform protein inference automatically to resolve protein ambiguity,
and some functions add a "gene" column by mapping protein accessions to gene symbols.
Besides, all `read_xxx()` functions now aggregate quantification from PSM to the glycopeptide level.

Detailed list:

* Add `read_pglyco3()` to read results from pGlyco3 (with built-in quantification).
* Add `read_byonic_byologic()` to read results from Byonic with Byologic quantification.
* All functions now automatically perform PSM aggregation to ensure glycopeptide-level quantification.
* Redesign `read_msfragger()` to support multiple PSM files. No awkward manul post-processing is needed anymore.

## Minor improvements

* Remove `protein_inference_method` parameter from `read_pglyco3_pglycoquant()`.
  The function now performs the "parsimony" protein inference algorithm.
* Optimize the logic of protein inference. For ties, alphabetical order is used.

# glyread 0.3.2

## Minor improvements

* Use `AnnotationDbi::mapIds()` to convert protein IDs to gene symbols.
* Remove `clusterProfiler` from Suggests.

## Bug fixes

* Try to fix the R CMD check errors on Github Actions.

# glyread 0.3.1

## Minor improvements

* `read_byonic_pglycoquant()` now only perform gene symbol conversion if both `clusterProfiler` 
   and the specified OrgDb package (e.g., `org.Hs.eg.db`) are installed. 
   If not, the gene column is omitted and a message is shown. 
   This improves robustness in minimal or CI environments.
* `clusterProfiler` and `org.Hs.eg.db` are now listed under Suggests, not Imports.

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
