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

test_that("read_glycan_finder extracts gene from protein accession", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N")

  # Gene should be extracted (format: P09871|C1S_HUMAN -> C1S_HUMAN)
  expect_true("gene" %in% names(exp$var_info))
  expect_true(all(stringr::str_detect(exp$var_info$gene, "^[A-Z0-9_]+$")))
})

test_that("read_glycan_finder extracts protein correctly", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N")

  # Protein should not contain the pipe symbol (isoform info removed)
  expect_false(any(stringr::str_detect(exp$var_info$protein, "[|]")))
})

test_that("read_glycan_finder extracts peptide_site correctly", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N")

  # peptide_site should be integer or numeric
  expect_true(is.integer(exp$var_info$peptide_site) || is.numeric(exp$var_info$peptide_site))
  # All values should be positive (site positions)
  expect_true(all(exp$var_info$peptide_site > 0, na.rm = TRUE))
})

test_that("read_glycan_finder filters by glycan_type correctly", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp_n <- read_glycan_finder(fp, glycan_type = "N")
  exp_o <- read_glycan_finder(fp, glycan_type = "O-GalNAc")

  # N-type should have different results from O-type
  expect_true(nrow(exp_n$var_info) > 0)
  # O-type should have fewer or equal rows
  expect_lte(nrow(exp_o$var_info), nrow(exp_n$var_info))
})
