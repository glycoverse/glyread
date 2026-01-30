test_that("read_glycan_finder() returns a glyexp experiment object", {
  result <- read_glycan_finder("data/glycan-finder-result.csv", glycan_type = "N")
  expect_s3_class(result, "glyexp_experiment")
})

test_that(".parse_glycan_finder_peptide() removes modifications", {
  peptide <- "NC(+57.02)GVN(+1913.68)C(+57.02)SGDVF"
  result <- glyread:::.parse_glycan_finder_peptide(peptide)
  expect_equal(result, "NCGVNCSGDVF")
})

test_that(".extract_glycan_finder_peptide_site() extracts N-glycan site", {
  peptide <- "NC(+57.02)GVN(+1913.68)C(+57.02)SGDVF"
  result <- glyread:::.extract_glycan_finder_peptide_site(peptide, "N")
  expect_equal(result, 5L)  # N at position 5 has the glycan modification
})

test_that(".extract_glycan_finder_peptide_site() extracts O-glycan site", {
  peptide <- "LGN(+2786.96)WSAMPS(+755.30)C(+57.02)K"
  result <- glyread:::.extract_glycan_finder_peptide_site(peptide, "O-GalNAc")
  expect_equal(result, 9L)  # S at position 9 has the glycan modification
})
