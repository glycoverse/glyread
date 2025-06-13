# ----- Label-free -----
test_that("it returns correct information (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      test_path("pglyco3-LFQ-sample-info.csv"),
      quant_method = "label-free"
    )
  )

  expect_equal(
    colnames(res$var_info),
    c(
      "variable", "peptide", "proteins", "genes", "glycan_composition",
      "glycan_structure", "peptide_site", "protein_sites", "charge",
      "modifications"
    )
  )
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  expect_s3_class(res$var_info$glycan_structure, "glyrepr_structure")
  expect_equal(colnames(res$sample_info), c("sample", "group"))
  expect_equal(colnames(res$expr_mat), c("S1", "S2"))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(res$meta_data, list(
    experiment_type = "glycoproteomics",
    glycan_type = "N",
    quantification_method = "label-free"
  ))
})


test_that("it accepts a sample_info tibble", {
  suppressMessages({
    sample_info <- readr::read_csv(test_path("pglyco3-LFQ-sample-info.csv"))
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      sample_info = sample_info,
      quant_method = "label-free"
    )
  })
  expect_equal(colnames(res$sample_info), c("sample", "group"))
})


test_that("it provides a default sample information tibble (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  expect_equal(res$sample_info, tibble::tibble(sample = c("S1", "S2")))
})


test_that("it raiss an error when samples in sample information is inconsistent with pGlyco3 result (label-free)", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  new_sample_info_path <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(sample_info, new_sample_info_path)
  expect_snapshot(
    read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      test_path(new_sample_info_path),
      quant_method = "label-free"
    ),
    error = TRUE
  )
})


test_that("zeros are replaced by NA (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  expect_true(sum(is.na(res$expr_mat)) > 0)
})


test_that("sample name converter works", {
  sample_name_converter <- function(x) {
    stringr::str_replace(x, "S", "Sample_")
  }
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      sample_name_converter = sample_name_converter
    )
  )
  expect_equal(colnames(res$expr_mat), c("Sample_1", "Sample_2"))
  expect_equal(res$sample_info$sample, c("Sample_1", "Sample_2"))
})
