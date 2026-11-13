# Original GlycReSoft search results

These are complete, unmodified rows from public **original identification CSV
exports**, not simulated search results or handwritten format examples.
Only selected groups were extracted; every member of a selected group is kept.
The R tests explicitly distinguish these source fixtures from derived malformed
or boundary-condition inputs created temporarily during a test.

Authors: Joshua Klein, Luis Carvalho, and Joseph Zaia.
Paper: *Expanding N-glycopeptide identifications by modeling fragmentation,
elution, and glycome connectivity*, Nature Communications (2024).
https://doi.org/10.1038/s41467-024-50338-5
Table 1 reports GlycReSoft 0.4.23; CSV files themselves contain no version column.
Do not infer exact exporter builds from publication dates.

Archive: https://doi.org/10.6084/m9.figshare.24578857.v1 (published 2024-07-23).
Metadata: https://api.figshare.com/v2/articles/24578857
License: **Creative Commons Attribution 4.0 International**, separate from
this package's MIT license: https://creativecommons.org/licenses/by/4.0/.
Attribution above and the description of row selection identify the source and
changes. No endorsement by the original authors is implied.

## Downloads and selections

- `PXD009654_dataset_base_glycopeptides.zip`:
  https://ndownloader.figshare.com/files/47117437
  MD5 `8f0fca7aebdac30f63af12c205c3210c`.
  - `20160301_Serum_SA_glycopeptide_HCC_2.glycopeptides.csv` -> `HCC_2.csv`:
    groups 1, 81, 118, 371.
  - `20160301_Serum_SA_glycopeptide_HCC_3.glycopeptides.csv` -> `HCC_3.csv`:
    groups 1, 289.
  - `20160301_Serum_SA_glycopeptide_Normal_3.glycopeptides.csv` -> `Normal_3.csv`:
    groups 1, 331, 372.
- `yeast_entrapment_dataset_base_glycopeptides.zip`:
  https://ndownloader.figshare.com/files/47116882
  MD5 `5b0790a2dc71eef22c22b649e4ff2731`.
  - `cwq_mix2-1_726.glycopeptides.csv` -> `cwq_mix2-1.csv`:
    groups 1, 481, 832 (shared proteins, Neu5Gc, and zero signal).

`inst/extdata/glycresoft.csv` is an exact copy of `HCC_2.csv` for runnable
examples; its provenance and license are also installed beside it.

## SHA-256

| File | Full original CSV | Fixture |
|---|---|---|
| HCC_2 | b05775d45638d660fe56e0778b72fd230c2e20c614488c84da4ab948afd49057 | f5f85b8657a6803b7ea6dd3cc0bc6ae0bbdd6c1782227e814f4d8723d0dd8198 |
| HCC_3 | f92e05c4f66f04bcda43b5ca97f55a872e574a14e2fdfdf4fa01d59c8696d101 | eacb020dd5270471e257339cfea92f317d030dd71cb112a0e8363d3e166bab75 |
| Normal_3 | 463ad8b22c95f77b2314fa99f5884bc6c87643f5365b86d7fdafa66273cb5ec2 | aa252a181a122fe59c5d52d08180dc3a7f38ba046f764e60a8e8bbefe548a5e4 |
| cwq_mix2-1 | 326df31708670c5037139549d8c7964f746589a2731817a06da970b4813038d2 | aae91784b15e9d0b587168b1dfda6b98add2362476a6f46ca6ee20976e0c98e9 |

## Semantics independently checked against official code

`group_id` is the database `ambiguous_id`. In `serialize/serializer.py`,
`save_glycopeptide_identification_set()` groups identification alternatives by
chromatogram before assigning an ambiguity ID. It is not a globally stable
feature ID. `output/csv_format.py` exports its value and `total_signal`.
The official 0.4.24 release has the same relevant schema as these files.
The master development schema is different and is not claimed as supported.
https://github.com/mobiusklein/glycresoft/commit/bd2ca5fd67b66ce0394cb2a7bef92f76cfa041ac
introduced `group_id` (included in v0.4.22).

Source peptide starts and site annotations are zero-based. Occupied sites are
read from the actual glycosylation markers, NOT every potential sequon listed
in `n_glycosylation_sites`. The raw localization probabilities remain available.

No old DB conversion, schema migration, manually invented CSV, or SQLite reader
is part of these fixtures or this implementation.
