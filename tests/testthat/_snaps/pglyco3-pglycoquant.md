# error when samples in sample information is inconsistent with pGlyco3 result (label-free)

    Code
      read_pglyco3_pglycoquant(test_path("pglyco3-pglycoquant-LFQ-result.list"),
      test_path(new_sample_info_path), name = "my_exp", quant_method = "label-free")
    Condition
      Error in `glyexp::experiment()`:
      ! Samples or variables must be consistent between `expr_mat`, `sample_info`, and `var_info`.
      x  Samples in `sample_info` but not in `expr_mat`: "S3"

# missing samples in sample_info file (TMT)

    Code
      read_pglyco3_pglycoquant(test_path("pglyco3-pglycoquant-TMT-result.list"),
      test_path(new_sample_info_path), name = "my_exp", quant_method = "TMT",
      tmt_type = "TMT-10plex", ref_channel = "126")
    Condition
      Error in `.read_pglyco3_pglycoquant_tmt()`:
      ! Sample information is missing for these pairs: "RAW2: 128N"

# extra samples in sample_info file (TMT)

    Code
      read_pglyco3_pglycoquant(test_path("pglyco3-pglycoquant-TMT-result.list"),
      test_path(new_sample_info_path), name = "my_exp", quant_method = "TMT",
      tmt_type = "TMT-10plex", ref_channel = "126")
    Condition
      Error in `.read_pglyco3_pglycoquant_tmt()`:
      ! Extra samples found in sample information: "S7", "S8", and "S9"

# intensities are normalized to reference channel

    Code
      res$expr_mat
    Output
                        S1  S2  S3  S4  S5  S6
      JNSDISSTR_1_H4N2 1.0 1.0 1.0 1.0 1.0 1.0
      JNSDISSTR_1_H6N2 1.5 1.5 1.5 1.5 1.5 1.5
      JNSDISSTR_1_H5N2 2.0 2.0 2.0 2.0 2.0 2.0

# duplicated samples in sample_info raises an error

    Code
      read_pglyco3_pglycoquant(test_path("pglyco3-pglycoquant-TMT-result.list"),
      test_path(new_sample_info_path), name = "my_exp", quant_method = "TMT",
      tmt_type = "TMT-10plex", ref_channel = "126")
    Condition
      Error in `.read_pglyco3_pglycoquant_tmt()`:
      ! Sample names in the sample information file must be unique.

