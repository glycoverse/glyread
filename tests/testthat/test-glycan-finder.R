test_that("read_glycan_finder reads basic CSV", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp <- read_glycan_finder(fp, glycan_type = "N")

  expect_s3_class(exp, "experiment")
  expect_equal(exp$glycan_type, "N")
})

test_that("read_glycan_finder filters by glycan_type", {
  fp <- test_path("data/glycan-finder-result.csv")
  exp_n <- read_glycan_finder(fp, glycan_type = "N")
  exp_o <- read_glycan_finder(fp, glycan_type = "O-GalNAc")

  # O-type should have fewer or equal rows
  expect_lte(nrow(exp_o$var_info), nrow(exp_n$var_info))
})
