# Original O-Pair output

`oglyco.psmtsv` is the unchanged MetaMorpheus project test fixture, containing
nine mouse-islet experimental PSM records, downloaded on 2026-11-13.

- Source: https://github.com/smith-chem-wisc/MetaMorpheus/blob/b8929b1ecde1ffc75f0fba4a61556eb4b155598e/MetaMorpheus/Test/TestData/oglyco.psmtsv
- Raw URL: https://raw.githubusercontent.com/smith-chem-wisc/MetaMorpheus/b8929b1ecde1ffc75f0fba4a61556eb4b155598e/MetaMorpheus/Test/TestData/oglyco.psmtsv
- SHA256: 78be7d5f3c3c9f57590f27d661c61e3f64c233a8226c578ccfad7db681ce93fb
- License: MIT; the upstream notice is included in LICENSE.txt.
- Format documentation: https://github.com/smith-chem-wisc/MetaMorpheus/wiki/O-Glyco-Search-Task
- Composition mapping: MetaMorpheus/EngineLayer/GlycoSearch/Glycan.cs (H/N/A/G/F).

This is an identification output, not a quantitative export. Matched-ion
intensities are fragment-spectrum measurements and must not become an abundance
assay. Reader regression tests start from these actual PSMs; altered copies
exercise filters and invalid-input cases, not invented exporter examples.
