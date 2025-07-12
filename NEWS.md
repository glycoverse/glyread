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
