# GproDIA original-result fixtures

These fixtures are row subsets of original GproDIA results deposited by the
GproDIA authors, not simulated exports.

Yang, Y., Yan, G., Kong, S., Wu, M., Yang, P., Cao, W., and Qiao, L. (2021).
GproDIA enables data-independent acquisition glycoproteomics with comprehensive
statistical control. Nature Communications 12, 6073.
https://doi.org/10.1038/s41467-021-26246-3

Archive: PXD023980 / IPX0002792000 / IPX0002792001
https://proteomecentral.proteomexchange.org/ui?pxid=PXD023980

Downloaded 2026-11-13 from:

- http://download.iprox.org/IPX0002792000/IPX0002792001/fissionyeast_aligned.tsv
- http://download.iprox.org/IPX0002792000/IPX0002792001/serum_uis_aligned.tsv

Original-file SHA-256:

- fissionyeast: `61175b2864e18df598ee922bfa4f64a3bec4e899a87e4a9534a7ce82b7f0e148`
- serum: `6fddf507944ec528413a8ede45d1be046fbe8cf7b273137b6c39241a8e92ae8e`

Changes: retain original header, all columns and all rows for selected
`transition_group_id` values; normalize line endings to LF. No cell values were
edited. Fission yeast IDs: 0, 1, 70, 75, 270, 1009. Serum IDs: 20, 32, 33,
410, 418, 498, 960, 995. These cover multiple runs, decoys, non-rank-1 aligned
features, floating-point negative q-values, multiple protein mappings, and
modified/unmodified precursors.

License: CC BY 4.0, https://creativecommons.org/licenses/by/4.0/ .
The iProX data policy explicitly assigns CC BY 4.0 to data submitted before
2023-03-23 (this deposit is from 2021):
https://www.iprox.cn/page/iproxDataLisence.html
This data license is separate from GproDIA's BSD-3-Clause software license.

Official descriptions of these files and export semantics:
https://github.com/lmsac/GproDIA/blob/1b11543de614fbb894476363f24f1385fc240134/docs/getting_started.md
https://github.com/lmsac/GproDIA/blob/1b11543de614fbb894476363f24f1385fc240134/docs/glycoform_inference.md
https://github.com/lmsac/GproDIA/blob/1b11543de614fbb894476363f24f1385fc240134/src/glycoprophet/export.py

Tests may deliberately mutate copies of these real rows to exercise validation;
those copies are not presented as original outputs.
