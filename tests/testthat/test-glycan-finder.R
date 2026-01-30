test_that("read_glycan_finder() returns a glyexp experiment object", {
  result <- read_glycan_finder("data/glycan-finder-result.csv", glycan_type = "N")
  expect_s3_class(result, "glyexp_experiment")
})

test_that(".parse_glycan_finder_peptide() removes modifications", {
  peptide <- "NC(+57.02)GVN(+1913.68)C(+57.02)SGDVF"
  result <- glyread:::.parse_glycan_finder_peptide(peptide)
  expect_equal(result, "NCGVNCSGDVF")
})
