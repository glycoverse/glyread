test_that("read_pglyco3 returns an experiment with correct information", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      test_path("pglyco3-sample-info.csv"),
      name = "my_exp"
    )
  )

  expect_equal(
    colnames(res$var_info),
    c("variable",
      "charge",
      "peptide",
      "modifications",
      "glycan_composition",
      "glycan_structure",
      "n_hex",
      "n_hexnac",
      "n_neuac",
      "n_fuc",
      "glycan_type",
      "bisecting",
      "n_antennae",
      "n_core_fuc",
      "n_arm_fuc",
      "n_gal",
      "n_terminal_gal",
      "peptide_site",
      "proteins",
      "genes",
      "protein_sites")
  )
  expect_equal(colnames(res$sample_info), c("sample", "group"))
  expect_equal(colnames(res$expr_mat), c("S1", "S2"))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(res$name, "my_exp")
  expect_true(!is.null(res$glycan_graphs))
})


test_that("read_pglyco3 provides a default name with current time", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      test_path("pglyco3-sample-info.csv"),
      parse_structure = FALSE
    )
  )

  expect_true(stringr::str_detect(
    res$name, "^exp \\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"
  ))
})


test_that("read_pglyco3 provides a default sample information tibble", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      name = "my_exp",
      parse_structure = FALSE
    )
  )
  expect_equal(res$sample_info, tibble::tibble(sample = c("S1", "S2")))
})


test_that("describe_glycans is always FALSE when parse_structure is FALSE", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      name = "my_exp",
      parse_structure = FALSE
    )
  )
  expect_true(!"glycan_type" %in% colnames(res$var_info))
})


test_that("quantification column depends on quantify_on argument", {
  suppressMessages(
    res1 <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      name = "my_exp",
      parse_structure = FALSE,
      quantify_on = "mono"
    )
  )
  expect_equal(res1$expr_mat[1, 1], 1e7)

  suppressMessages(
    res2 <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      name = "my_exp",
      parse_structure = FALSE,
      quantify_on = "sum"
    )
  )
  expect_equal(res2$expr_mat[1, 1], 1.5e7)
})


test_that("duplicated records are summed", {
  suppressMessages(df <- readr::read_delim(test_path("pglyco3-result.txt")))
  df <- dplyr::bind_rows(df, df)
  res <- withr::with_tempdir({
    readr::write_delim(df, "pglyco3-result.txt", delim = "\t")
    suppressMessages(
      read_pglyco3(
        "pglyco3-result.txt",
        name = "my_exp",
        parse_structure = FALSE
      )
    )
  })

  expect_snapshot(res$expr_mat)
})


test_that("glycan compositions are correct", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      test_path("pglyco3-sample-info.csv"),
      name = "my_exp"
    )
  )

  expected <- tibble::tribble(
    ~n_hex, ~n_hexnac, ~n_neuac, ~n_fuc,
    4, 4, 0, 1,
    5, 2, 0, 0,
    3, 2, 0, 0,
    5, 4, 0, 1,
    5, 4, 1, 1,
    3, 2, 0, 0
  )
  expect_equal(res$var_info[c("n_hex", "n_hexnac", "n_neuac", "n_fuc")], expected)
})


test_that("differ_a_g setting to FALSE changes A to S", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("pglyco3-result.txt"),
      name = "my_exp",
      differ_a_g = FALSE
    )
  )

  expected <- c("H4N4F1", "H5N2", "H3N2", "H5N4F1", "H5N4F1S1", "H3N2")
  expect_equal(res$var_info$glycan_composition, expected)
  expect_snapshot(res$glycan_graphs[[5]])
})


test_that("error when samples in sample information is inconsistent with pGlyco3 result", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  withr::with_tempdir({
    readr::write_csv(sample_info, "sample_info.csv")
    expect_error(
      read_pglyco3(
        test_path("sample_info.csv"),
        name = "my_exp",
        differ_a_g = FALSE
      )
    )
  })
})
