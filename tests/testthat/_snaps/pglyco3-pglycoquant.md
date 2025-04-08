# it raiss an error when samples in sample information is inconsistent with pGlyco3 result (label-free)

    Code
      read_pglyco3_pglycoquant(test_path("pglyco3-pglycoquant-LFQ-result.list"),
      test_path(new_sample_info_path), name = "my_exp", quant_method = "label-free")
    Condition
      Error in `glyexp::experiment()`:
      ! Samples or variables must be consistent between `expr_mat`, `sample_info`, and `var_info`.
      x  Samples in `sample_info` but not in `expr_mat`: "S3"

