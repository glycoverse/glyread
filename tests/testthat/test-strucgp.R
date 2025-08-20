test_that("read_strucgp works", {
  result <- read_strucgp(test_path("data/strucgp-result.xlsx"))
  
  # Check that result is a tibble
  expect_s3_class(result, "tbl_df")
  
  # Check expected columns are present
  expected_cols <- c(
    "peptide", "protein", "gene", "protein_site", 
    "glycan_composition", "glycan_structure"
  )
  expect_true(all(expected_cols %in% colnames(result)))
  
  # Check column types
  expect_type(result$peptide, "character")
  expect_type(result$protein, "character") 
  expect_type(result$gene, "character")
  expect_type(result$protein_site, "character")
  expect_s3_class(result$glycan_composition, "glyrepr_composition")
  expect_s3_class(result$glycan_structure, "glyrepr_structure")
  
  # Check that result has data
  expect_gt(nrow(result), 0)
  
  # Check that duplicates are removed (distinct() is applied)
  expect_equal(nrow(result), nrow(dplyr::distinct(result)))
})

test_that(".parse_strucgp_comp works correctly", {
  # Test the internal composition parsing function
  test_comp <- "H5N4F1S2+Na"
  
  # Access the internal function 
  result <- glyread:::.parse_strucgp_comp(test_comp)
  
  expect_s3_class(result, "glyrepr_composition")
  
  # Check that the composition string is parsed correctly
  # H -> Hex, N -> HexNAc, F -> dHex, S -> NeuAc
  # +Na should be removed
  expected_comp <- glyrepr::as_glycan_composition("Hex(5)HexNAc(4)dHex(1)NeuAc(2)")
  expect_equal(result, expected_comp)
})

test_that(".parse_strucgp_comp handles edge cases", {
  # Test with no modifications
  test_comp1 <- "H3N2"
  result1 <- glyread:::.parse_strucgp_comp(test_comp1)
  expected1 <- glyrepr::as_glycan_composition("Hex(3)HexNAc(2)")
  expect_equal(result1, expected1)
  
  # Test with single digit numbers
  test_comp2 <- "H1N1F1S1"
  result2 <- glyread:::.parse_strucgp_comp(test_comp2)
  expected2 <- glyrepr::as_glycan_composition("Hex(1)HexNAc(1)dHex(1)NeuAc(1)")
  expect_equal(result2, expected2)
})
