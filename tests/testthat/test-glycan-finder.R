test_that("read_glycan_finder() returns a glyexp experiment object", {
  result <- read_glycan_finder("data/glycan-finder-result.csv", glycan_type = "N")
  expect_s3_class(result, "glyexp_experiment")
})
