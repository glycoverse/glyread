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
  result <- glyread:::.tidy_glycan_finder(df, "N")

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
