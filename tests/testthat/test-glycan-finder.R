test_that("read_glycan_finder() returns a glyexp experiment object", {
  result <- suppressMessages(read_glycan_finder("data/glycan-finder-result.csv", glycan_type = "N"))
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

test_that(".parse_glycan_finder_protein() extracts first accession", {
  protein <- "P09871|C1S_HUMAN"
  result <- glyread:::.parse_glycan_finder_protein(protein)
  expect_equal(result, "P09871")
})

test_that(".parse_glycan_finder_composition() parses composition string", {
  comp <- "(HexNAc)4(Hex)5(NeuAc)1"
  result <- glyread:::.parse_glycan_finder_composition(comp)
  expect_s3_class(result, "glyrepr_composition")
  expect_equal(unname(unlist(result)), c(5L, 4L, 1L))
})

test_that(".read_glycan_finder_df() reads CSV correctly", {
  result <- glyread:::.read_glycan_finder_df("data/glycan-finder-result.csv")
  expect_s3_class(result, "tbl_df")
  expect_true("Protein Accession" %in% colnames(result))
  expect_true("Peptide" %in% colnames(result))
  expect_true("Glycan" %in% colnames(result))
  expect_true("Glycan Type" %in% colnames(result))
})

test_that(".filter_glycan_finder_by_type() filters correctly", {
  df <- tibble::tribble(
    ~`Glycan Type`, ~value,
    "N-Link", 1,
    "O-Link", 2,
    "N-Link;O-Link", 3
  )

  result_n <- glyread:::.filter_glycan_finder_by_type(df, "N")
  expect_equal(nrow(result_n), 2)  # N-Link and N-Link;O-Link

  result_o <- glyread:::.filter_glycan_finder_by_type(df, "O-GalNAc")
  expect_equal(nrow(result_o), 2)  # O-Link and N-Link;O-Link
})

test_that(".tidy_glycan_finder() transforms data correctly", {
  df <- glyread:::.read_glycan_finder_df("data/glycan-finder-result.csv")
  result <- suppressMessages(glyread:::.tidy_glycan_finder(df, "N"))

  expect_s3_class(result, "tbl_df")
  expect_true("peptide" %in% colnames(result))
  expect_true("protein" %in% colnames(result))
  expect_true("peptide_site" %in% colnames(result))
  expect_true("protein_site" %in% colnames(result))
  expect_true("glycan_composition" %in% colnames(result))
  expect_true("glycan_structure" %in% colnames(result))
  expect_true("sample" %in% colnames(result))
  expect_true("value" %in% colnames(result))
})

test_that(".select_glycan_element() selects correct element for mixed types", {
  # Test N-glycan selection from mixed type
  glycan_types <- c("N-Link;O-Link", "N-Link;O-Link", "N-Link")
  glycans <- c("(HexNAc)4(Hex)5(NeuAc)4;(HexNAc)3(Fuc)1", "(HexNAc)7(Hex)4(NeuAc)1;(HexNAc)1", "(HexNAc)4(Hex)5(NeuAc)1")

  result_n <- glyread:::.select_glycan_element(glycan_types, glycans, "N")
  expect_equal(result_n[1], "(HexNAc)4(Hex)5(NeuAc)4")  # First element for mixed
  expect_equal(result_n[2], "(HexNAc)7(Hex)4(NeuAc)1")  # First element for mixed
  expect_equal(result_n[3], "(HexNAc)4(Hex)5(NeuAc)1")  # Unchanged for single type

  result_o <- glyread:::.select_glycan_element(glycan_types, glycans, "O-GalNAc")
  expect_equal(result_o[1], "(HexNAc)3(Fuc)1")  # Second element for mixed
  expect_equal(result_o[2], "(HexNAc)1")        # Second element for mixed
  expect_equal(result_o[3], "(HexNAc)4(Hex)5(NeuAc)1")  # Unchanged for single type
})

# Regression tests for previous fixes

test_that("PTM-based detection correctly identifies small glycans like O-GlcNAc", {
  # O-GlcNAc is ~203 Da, which would be missed by a mass > 500 threshold
  peptide <- "N(+2360.86)ATSSSSQDPES(+203.08)LQDR"
  ptm <- "(HexNAc)7(Hex)4(NeuAc)1; (HexNAc)1"

  # The O-glycan (O-GlcNAc) is on the S with +203.08 (position 12)
  result <- glyread:::.extract_glycan_finder_peptide_site(peptide, "O-GalNAc", ptm)
  expect_equal(result, 12L)  # S at position 12

  # The N-glycan is on the N with +2360.86 (position 1)
  result_n <- glyread:::.extract_glycan_finder_peptide_site(peptide, "N", ptm)
  expect_equal(result_n, 1L)  # N at position 1
})

test_that("sample columns exclude summary columns like Area C3 and Area H", {
  # Use full pipeline to get processed data
  result <- suppressMessages(read_glycan_finder("data/glycan-finder-result.csv", glycan_type = "N"))

  # Get sample info and check columns
  sample_info <- glyexp::get_sample_info(result)
  sample_names <- sample_info$sample

  # Should only have replicate columns (C_3, H_3, H_1, H_2), not summary columns (C3, H)
  expect_true("C_3" %in% sample_names)
  expect_true("H_3" %in% sample_names)
  expect_true("H_1" %in% sample_names)
  expect_true("H_2" %in% sample_names)
  expect_false("C3" %in% sample_names)  # Summary column should be excluded
  expect_false("H" %in% sample_names)   # Summary column should be excluded
})

test_that("end-to-end with mixed glycan types does not confuse N and O glycans", {
  # Read the test data - it contains mixed type rows like:
  # "N-Link;O-Link" with Glycan "(HexNAc)4(Hex)5(NeuAc)4;(HexNAc)3(Fuc)1"
  result_n <- suppressMessages(read_glycan_finder("data/glycan-finder-result.csv", glycan_type = "N"))
  result_o <- suppressMessages(read_glycan_finder("data/glycan-finder-result.csv", glycan_type = "O-GalNAc"))

  var_info_n <- glyexp::get_var_info(result_n)
  var_info_o <- glyexp::get_var_info(result_o)

  # When reading N-glycans, we should NOT see pure O-glycan compositions
  # (e.g., single HexNAc without Hex)
  n_comps <- var_info_n$glycan_composition
  o_comps <- var_info_o$glycan_composition

  # Check that N-glycans and O-glycans are different sets
  n_comp_str <- purrr::map_chr(n_comps, ~ as.character(.x))
  o_comp_str <- purrr::map_chr(o_comps, ~ as.character(.x))

  # There should be some compositions unique to each type
  expect_false(length(setdiff(n_comp_str, o_comp_str)) == 0)
  expect_false(length(setdiff(o_comp_str, n_comp_str)) == 0)
})
