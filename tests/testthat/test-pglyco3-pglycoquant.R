test_that("read_pglyco3_pglycoquant returns an experiment with correct information", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-result.list"),
      test_path("pglyco3-sample-info.csv"),
      name = "my_exp"
    )
  )

  expect_equal(
    colnames(res$var_info),
    c("variable",
      "peptide",
      "proteins",
      "genes",
      "glycan_composition",
      "glycan_structure",
      "peptide_site",
      "protein_sites",
      "charge",
      "modifications")
  )
  expect_equal(colnames(res$sample_info), c("sample", "group"))
  expect_equal(colnames(res$expr_mat), c("S1", "S2"))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(res$name, "my_exp")
  expect_true(is.null(res$glycan_graphs))
})


test_that("read_pglyco3_pglycoquant returns an experiment with correct information with structures parsed", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-result.list"),
      test_path("pglyco3-sample-info.csv"),
      name = "my_exp",
      parse_structure = TRUE,
    )
  )

  expect_true(!is.null(res$glycan_graphs))
})


test_that("read_pglyco3_pglycoquant provides a default name with current time", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-result.list"),
      test_path("pglyco3-sample-info.csv"),
    )
  )

  expect_true(stringr::str_detect(
    res$name, "^exp \\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"
  ))
})


test_that("read_pglyco3_pglycoquant provides a default sample information tibble", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-result.list"),
      name = "my_exp",
    )
  )
  expect_equal(res$sample_info, tibble::tibble(sample = c("S1", "S2")))
})


test_that("describe_glycans is always FALSE when parse_structure is FALSE", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-result.list"),
      name = "my_exp",
    )
  )
  expect_true(!"glycan_type" %in% colnames(res$var_info))
})


test_that("error when samples in sample information is inconsistent with pGlyco3 result", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  withr::with_tempdir({
    readr::write_csv(sample_info, "sample_info.csv")
    expect_error(
      read_pglyco3_pglycoquant(
        test_path("pglyco3-pglycoquant-result.list"),
        test_path("sample_info.csv"),
        name = "my_exp",
      )
    )
  })
})


test_that("zeros are replaced by NA", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-result.list"),
      name = "my_exp",
    )
  )
  expect_true(sum(is.na(res$expr_mat)) > 0)
})
