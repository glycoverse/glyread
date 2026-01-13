test_that("read_strucgp works", {
  suppressMessages(result <- read_strucgp(test_path("data/strucgp-result.xlsx")))

  # Check that result is an experiment object
  expect_s3_class(result, "glyexp_experiment")

  # Check var_info has expected columns (peptide is removed after standardization)
  expected_cols <- c(
    "variable", "protein", "gene", "protein_site",
    "glycan_composition", "glycan_structure"
  )
  expect_true(all(expected_cols %in% colnames(result$var_info)))

  # Check column types in var_info
  expect_type(result$var_info$variable, "character")
  expect_type(result$var_info$protein, "character")
  expect_type(result$var_info$gene, "character")
  expect_type(result$var_info$protein_site, "integer")
  expect_s3_class(result$var_info$glycan_composition, "glyrepr_composition")
  expect_s3_class(result$var_info$glycan_structure, "glyrepr_structure")

  # Check that var_info has data
  expect_gt(nrow(result$var_info), 0)

  # Check that expr_mat is a binary (0/1) matrix
  expect_true(all(result$expr_mat %in% c(0, 1)))
  expect_equal(nrow(result$expr_mat), nrow(result$var_info))

  # Check that at least some glycopeptides were identified (have 1s)
  expect_true(any(result$expr_mat == 1))

  # Check that sample_info exists
  expect_s3_class(result$sample_info, "tbl_df")
  expect_true("sample" %in% colnames(result$sample_info))
  expect_gt(nrow(result$sample_info), 0)

  # Check dimensions match
  expect_equal(ncol(result$expr_mat), nrow(result$sample_info))
})

test_that("read_strucgp works with parse_structure = FALSE", {
  suppressMessages(result <- read_strucgp(test_path("data/strucgp-result.xlsx"), parse_structure = FALSE))
  
  # Check that result is an experiment object
  expect_s3_class(result, "glyexp_experiment")
  
  # Check var_info does NOT have glycan_structure column
  expect_false("glycan_structure" %in% colnames(result$var_info))
  
  # Check that glycan_composition is still present
  expect_true("glycan_composition" %in% colnames(result$var_info))
  expect_s3_class(result$var_info$glycan_composition, "glyrepr_composition")
})

test_that("read_strucgp works with custom sample_info", {
  # Create a simple sample_info with the actual sample name from the test file
  sample_info <- data.frame(
    sample = c("20250403-LZE-NGlyco5_3"),
    group = c("treatment")
  )
  
  suppressMessages(result <- read_strucgp(
    test_path("data/strucgp-result.xlsx"),
    sample_info = sample_info
  ))
  
  # Check that sample_info is included
  expect_true("group" %in% colnames(result$sample_info))
  expect_equal(as.character(result$sample_info$group), c("treatment"))
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
