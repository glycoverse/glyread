# ----- Label-free -----
test_that("read_pglyco3_pglycoquant returns correct information (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      test_path("pglyco3-LFQ-sample-info.csv"),
      name = "my_exp",
      quant_method = "label-free"
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


test_that("read_pglyco3_pglycoquant parses structures", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      test_path("pglyco3-LFQ-sample-info.csv"),
      name = "my_exp",
      quant_method = "label-free",
      parse_structure = TRUE
    )
  )

  expect_true(!is.null(res$glycan_graphs))
})


test_that("read_pglyco3_pglycoquant sums duplicated records (label-free)", {
  df <- suppressMessages(readr::read_tsv(test_path("pglyco3-pglycoquant-LFQ-result.list")))
  df <- dplyr::bind_rows(df, df)
  withr::with_tempdir({
    readr::write_tsv(df, "pglyco3-pglycoquant-result2.list")
    suppressMessages(
      exp <- read_pglyco3_pglycoquant(
        test_path("pglyco3-pglycoquant-result2.list"),
        name = "my_exp",
        quant_method = "label-free"
    ))
  })

  expect_equal(nrow(exp$var_info), nrow(df) / 2)
})


test_that("read_pglyco3_pglycoquant provides a default name with current time", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      test_path("pglyco3-LFQ-sample-info.csv"),
      quant_method = "label-free"
    )
  )

  expect_true(stringr::str_detect(
    res$name, "^exp \\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"
  ))
})


test_that("read_pglyco3_pglycoquant provides a default sample information tibble (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      name = "my_exp",
      quant_method = "label-free"
    )
  )
  expect_equal(res$sample_info, tibble::tibble(sample = c("S1", "S2")))
})


test_that("error when samples in sample information is inconsistent with pGlyco3 result (label-free)", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  withr::with_tempdir({
    readr::write_csv(sample_info, "sample_info.csv")
    expect_error(
      read_pglyco3_pglycoquant(
        test_path("pglyco3-pglycoquant-LFQ-result.list"),
        test_path("sample_info.csv"),
        name = "my_exp",
        quant_method = "label-free"
      )
    )
  })
})


test_that("zeros are replaced by NA (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-LFQ-result.list"),
      name = "my_exp",
      quant_method = "label-free"
    )
  )
  expect_true(sum(is.na(res$expr_mat)) > 0)
})


# ----- TMT -----
test_that("read_pglyco3_pglycoquant returns correct information (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-TMT-result.list"),
      test_path("pglyco3-TMT-sample-info.csv"),
      name = "my_exp",
      quant_method = "TMT",
      tmt_type = "TMT-10plex",
      ref_channel = "126"
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
  expect_equal(colnames(res$sample_info), c("raw_name", "channel", "sample", "group"))
  expect_equal(colnames(res$expr_mat), paste0("S", 1:6))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(res$name, "my_exp")
  expect_true(is.null(res$glycan_graphs))
})


test_that("missing samples in sample_info file (TMT)", {
  sample_info <- suppressMessages(readr::read_csv(test_path("pglyco3-TMT-sample-info.csv")))
  sample_info <- dplyr::slice_head(sample_info, n = 5)  # remove the last sample
  withr::with_tempdir({
    readr::write_csv(sample_info, "sample_info.csv")
    expect_error(
      read_pglyco3_pglycoquant(
        test_path("pglyco3-pglycoquant-TMT-result.list"),
        test_path("sample_info.csv"),
        name = "my_exp",
        quant_method = "TMT",
        tmt_type = "TMT-10plex",
        ref_channel = "126"
      )
    )
  })
})


test_that("extra samples in sample_info file (TMT)", {
  sample_info <- suppressMessages(readr::read_csv(test_path("pglyco3-TMT-sample-info.csv")))
  new_samples_df <- tibble::tribble(
    ~raw_name, ~channel, ~sample, ~group,
    "RAW3", "127N", "S7", "A",
    "RAW3", "127C", "S8", "B",
    "RAW3", "128N", "S9", "A",
  )
  sample_info <- dplyr::bind_rows(sample_info, new_samples_df)
  withr::with_tempdir({
    readr::write_csv(sample_info, "sample_info.csv")
    expect_error(
      read_pglyco3_pglycoquant(
        test_path("pglyco3-pglycoquant-TMT-result.list"),
        test_path("sample_info.csv"),
        name = "my_exp",
        quant_method = "TMT",
        tmt_type = "TMT-10plex",
        ref_channel = "126"
      )
    )
  })
})


test_that("intensities are normalized to reference channel", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-TMT-result.list"),
      test_path("pglyco3-TMT-sample-info.csv"),
      name = "my_exp",
      quant_method = "TMT",
      tmt_type = "TMT-10plex",
      ref_channel = "126"
    )
  )

  expect_snapshot(res$expr_mat)
})


test_that("read_pglyco3_pglycoquant provides a default sample information tibble (TMT)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("pglyco3-pglycoquant-TMT-result.list"),
      name = "my_exp",
      quant_method = "TMT",
      tmt_type = "TMT-10plex",
      ref_channel = "126"
    )
  )
  expected <- tibble::tribble(
    ~raw_name, ~channel, ~sample,
    "RAW1", "127N", "S1",
    "RAW2", "127N", "S2",
    "RAW1", "127C", "S3",
    "RAW2", "127C", "S4",
    "RAW1", "128N", "S5",
    "RAW2", "128N", "S6",
  )
  expect_equal(res$sample_info, expected)
})


test_that("duplicated samples in sample_info raises an error", {
  sample_info <- tibble::tribble(
    ~raw_name, ~channel, ~sample,
    "RAW1", "127N", "S1",
    "RAW2", "127N", "S2",
    "RAW1", "127C", "S3",
    "RAW2", "127C", "S4",
    "RAW1", "128N", "S5",
    "RAW2", "128N", "S5",  # duplicated sample
  )
  withr::with_tempdir({
    readr::write_csv(sample_info, "sample_info.csv")
    expect_error(
      read_pglyco3_pglycoquant(
        test_path("pglyco3-pglycoquant-TMT-result.list"),
        test_path("sample_info.csv"),
        name = "my_exp",
        quant_method = "TMT",
        tmt_type = "TMT-10plex",
        ref_channel = "126"
      )
    )
  })
})
