test_that("all-missing abundance does not warn", {
  tidy_df <- tibble::tibble(
    sample = "S1",
    value = 0,
    protein = "P1",
    protein_site = 1L,
    glycan_composition = "Hex(1)"
  )

  expect_no_warning(
    result <- suppressMessages(glyread:::.read_template(
      tidy_df,
      sample_info_arg = NULL,
      glycan_type = "N",
      quant_method = "label-free",
      composition_parser = glyrepr::as_glycan_composition,
      parse_structure = FALSE
    ))
  )

  expect_s4_class(result, "GlycoproteomicSE")
  expect_true(all(is.na(SummarizedExperiment::assay(result, "abundance"))))
})
