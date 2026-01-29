test_that("read_glycan_finder reads basic CSV", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N")

  expect_s3_class(exp, "glyexp_experiment")
  expect_equal(exp$meta_data$glycan_type, "N")
})

test_that("read_glycan_finder filters by glycan_type", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp_n <- read_glycan_finder(fp, glycan_type = "N")
  exp_o <- read_glycan_finder(fp, glycan_type = "O-GalNAc")

  # O-type should have fewer or equal rows
  expect_lte(nrow(exp_o$var_info), nrow(exp_n$var_info))
})

test_that("read_glycan_finder parses structures", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N", parse_structure = TRUE)

  expect_true("glycan_structure" %in% names(exp$var_info))
  expect_s3_class(exp$var_info$glycan_structure, "glyrepr_structure")
})

test_that("read_glycan_finder skips structure parsing when disabled", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N", parse_structure = FALSE)

  expect_false("glycan_structure" %in% names(exp$var_info))
})

test_that("read_glycan_finder replaces 0 with NA", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N")

  # Check that expression matrix has no zeros (they should be NA)
  expect_true(all(is.na(exp$expr_mat) | exp$expr_mat > 0))
})

test_that("read_glycan_finder applies sample_name_converter", {
  fp <- test_path("data/glycan-finder-result.csv")

  # Simple converter: prefix with "sample_"
  exp <- read_glycan_finder(
    fp,
    glycan_type = "N",
    sample_name_converter = function(x) paste0("sample_", x)
  )

  expect_true(all(stringr::str_starts(colnames(exp$expr_mat), "sample_")))
})
