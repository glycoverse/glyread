# Test UniProt identifier parsing functions

# ----- Test byonic-pglycoquant UniProt parsing -----
test_that(".convert_byonic_columns correctly extracts UniProt accessions", {
  # Create test data with various UniProt formats
  test_df <- tibble::tibble(
    `Protein Name` = c(
      ">sp|P08185|CBG_HUMAN",           # Standard format
      ">sp|P08185-1|CBG_HUMAN",         # Isoform format
      ">tr|A0A024R6P0|A0A024R6P0_HUMAN", # TrEMBL format
      ">tr|Q9Y6R7-2|FCGBP_HUMAN",       # TrEMBL with isoform
      ">sp|P19652|A1AG2_HUMAN"          # Another standard format
    ),
    Peptide = c("K.nGTR.G", "K.nGTR.G", "K.nGTR.G", "K.nGTR.G", "K.nGTR.G"),
    Composition = rep("HexNAc(4)Hex(4)Fuc(1)NeuAc(1)", 5),
    Position = rep(123L, 5)
  )
  
  # Apply the conversion function
  result <- glyread:::.convert_byonic_columns(test_df)
  
  # Check that accessions are correctly extracted
  expected_proteins <- c("P08185", "P08185-1", "A0A024R6P0", "Q9Y6R7-2", "P19652")
  expect_equal(result$protein, expected_proteins)
  
  # Check that isoform suffixes are preserved
  expect_true("P08185-1" %in% result$protein)
  expect_true("Q9Y6R7-2" %in% result$protein)
})

# ----- Test pglyco3 UniProt parsing -----
test_that(".convert_pglyco3_columns correctly extracts UniProt accessions", {
  # Create test data with various UniProt formats
  test_df <- tibble::tibble(
    Peptide = rep("TESTPEPTIDE", 4),
    Proteins = c(
      "sp|P08185|CBG_HUMAN",           # Standard format
      "sp|P08185-1|CBG_HUMAN",         # Isoform format
      "sp|P19652|A1AG2_HUMAN;sp|Q9Y6R7-2|FCGBP_HUMAN", # Multiple proteins with isoform
      "sp|A0A024R6P0|A0A024R6P0_HUMAN" # Another standard format
    ),
    Genes = rep("TESTGENE", 4),
    GlycanComposition = rep("H(4)N(4)F(1)A(1)", 4),
    PlausibleStruct = rep("(N(F)(N(H(H(N))(H(N(H))))))", 4),
    GlySite = rep(1L, 4),
    ProSites = rep("123", 4)
  )
  
  # Apply the conversion function
  result <- glyread:::.convert_pglyco3_columns(test_df)
  
  # Check that accessions are correctly extracted
  expected_proteins <- c(
    "P08185",
    "P08185-1", 
    "P19652;Q9Y6R7-2",
    "A0A024R6P0"
  )
  expect_equal(result$proteins, expected_proteins)
  
  # Check that isoform suffixes are preserved
  expect_true(any(stringr::str_detect(result$proteins, "P08185-1")))
  expect_true(any(stringr::str_detect(result$proteins, "Q9Y6R7-2")))
})

# ----- Test byonic-byologic UniProt parsing -----
test_that(".refine_byonic_byologic_columns correctly extracts UniProt accessions", {
  # Create test data with various UniProt formats
  test_df <- tibble::tibble(
    protein_name = c(
      "sp|P08185|CBG_HUMAN",           # Standard format
      "sp|P08185-1|CBG_HUMAN",         # Isoform format
      "tr|A0A024R6P0|A0A024R6P0_HUMAN", # TrEMBL format
      "tr|Q9Y6R7-2|FCGBP_HUMAN"        # TrEMBL with isoform
    ),
    sequence = rep("K.TESTPEPTIDE.G", 4),
    glycans = rep("HexNAc(4)Hex(4)Fuc(1)NeuAc(1)", 4),
    mod_summary = rep("N123(NGlycan)", 4),
    start_aa = rep(100L, 4),
    ms_alias_name = rep("sample1", 4),
    xic_area_summed = rep(1000.0, 4)
  )
  
  # Apply the conversion function
  result <- glyread:::.refine_byonic_byologic_columns(test_df)
  
  # Check that accessions are correctly extracted
  expected_proteins <- c("P08185", "P08185-1", "A0A024R6P0", "Q9Y6R7-2")
  expect_equal(result$protein, expected_proteins)
  
  # Check that isoform suffixes are preserved
  expect_true("P08185-1" %in% result$protein)
  expect_true("Q9Y6R7-2" %in% result$protein)
})

# ----- Test edge cases -----
test_that("UniProt parsing handles edge cases correctly", {
  # Test empty and malformed inputs for byonic-pglycoquant
  test_df_byonic <- tibble::tibble(
    `Protein Name` = c(
      ">sp|P08185|CBG_HUMAN",
      ">malformed_entry",
      ">sp||EMPTY_ACCESSION",
      NA_character_
    ),
    Peptide = rep("K.nGTR.G", 4),
    Composition = rep("HexNAc(4)Hex(4)Fuc(1)NeuAc(1)", 4),
    Position = rep(123L, 4)
  )
  
  result_byonic <- glyread:::.convert_byonic_columns(test_df_byonic)
  
  # Should extract valid accession and handle malformed entries gracefully
  expect_equal(result_byonic$protein[1], "P08185")
  expect_true(is.na(result_byonic$protein[2]) || result_byonic$protein[2] == "")
  expect_true(is.na(result_byonic$protein[3]) || result_byonic$protein[3] == "")
  expect_true(is.na(result_byonic$protein[4]))
})
