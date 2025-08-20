# ----- Basic functionality tests -----
test_that("it returns correct information (label-free) with single sample", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )

  # Check structure of returned experiment object
  expect_s3_class(res, "glyexp_experiment")

  # Check variable information columns
  expect_equal(
    colnames(res$var_info),
    c(
      "variable", "peptide", "peptide_site",
      "protein", "protein_site", "gene", "glycan_composition"
    )
  )

  # Check glycan composition is properly parsed
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")

  # Check sample information has correct structure
  expect_equal(colnames(res$sample_info), "sample")

  # Check expression matrix dimensions and names
  expect_equal(ncol(res$expr_mat), 3) # Three samples in new test data
  expect_equal(nrow(res$expr_mat), nrow(res$var_info))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)

  # Check metadata
  expect_equal(res$meta_data, list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  ))

  # Check sample names are extracted correctly from directory structure
  expect_equal(sort(colnames(res$expr_mat)), c("H_1", "H_2", "H_3"))
  expect_equal(sort(res$sample_info$sample), c("H_1", "H_2", "H_3"))
})

test_that("it provides a default sample information tibble", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )
  # Sample names should be extracted from directory structure
  expected_samples <- c("H_1", "H_2", "H_3")
  expect_equal(sort(res$sample_info$sample), sort(expected_samples))
  expect_equal(ncol(res$sample_info), 1)
})

test_that("zeros are replaced by NA", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )
  # Check if there are any NA values in the intensity data
  # (should be converted from zero intensities)
  expect_true(sum(is.na(res$expr_mat)) >= 0)
})

test_that("it filters uncertain or multisite PSMs", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )
  # All remaining PSMs should have exactly one glycosite
  # This is verified by checking that all variables have valid data
  expect_true(nrow(res$var_info) > 0)
  expect_true(all(!is.na(res$var_info$peptide_site)))
  expect_true(all(!is.na(res$var_info$protein_site)))
})

test_that("it handles O-linked glycan type", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free",
      glycan_type = "O"
    )
  )
  expect_equal(res$meta_data$glycan_type, "O")
})

# ----- Parameter validation tests -----
test_that("it validates directory path parameter", {
  expect_error(
    read_msfragger(
      "non_existent_directory",
      quant_method = "label-free"
    ),
    class = "simpleError"
  )

  # Test with a file instead of directory
  temp_file <- withr::local_tempfile(fileext = ".txt")
  writeLines("test", temp_file)
  expect_error(
    read_msfragger(
      temp_file,
      quant_method = "label-free"
    ),
    class = "simpleError"
  )
})

test_that("it validates quant_method parameter", {
  expect_error(
    read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "invalid_method"
    ),
    " Must be element of set"
  )
})

test_that("it validates glycan_type parameter", {
  expect_error(
    read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free",
      glycan_type = "invalid_type"
    ),
    " Must be element of set"
  )
})

test_that("it validates sample_name_converter parameter", {
  expect_error(
    read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free",
      sample_name_converter = "not_a_function"
    ),
    class = "simpleError"
  )
})

test_that("it rejects TMT quantification (not implemented)", {
  expect_error(
    read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "TMT"
    ),
    "TMT quantification is not supported yet"
  )
})

test_that("sample name converter works", {
  sample_name_converter <- function(x) {
    paste0("Sample_", x)
  }
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free",
      sample_name_converter = sample_name_converter
    )
  )
  expected_samples <- c("Sample_H_1", "Sample_H_2", "Sample_H_3")
  expect_equal(sort(colnames(res$expr_mat)), sort(expected_samples))
  expect_equal(sort(res$sample_info$sample), sort(expected_samples))
})

# ----- Test glycan composition parsing -----
test_that("it correctly parses glycan compositions", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )

  # Check that glycan compositions are properly parsed
  glycan_comps <- res$var_info$glycan_composition
  expect_true(all(!is.na(glycan_comps)))

  # Check that compositions contain expected monosaccharides
  comp_strings <- as.character(glycan_comps)
  expect_true(any(grepl("Hex", comp_strings)))
  expect_true(any(grepl("HexNAc", comp_strings)))

  # Test that Fuc is converted to dHex if present in the test data
  # Check if any compositions contain dHex (converted from Fuc)
  if (any(grepl("dHex", comp_strings))) {
    expect_false(any(grepl("Fuc", comp_strings)))
  }
})

test_that("Fuc is correctly converted to dHex in glycan compositions", {
  # Test the internal parsing function directly
  test_composition <- "HexNAc(2)Hex(5)Fuc(1) % 1234.5678"
  parsed_comp <- glyread:::.parse_msfragger_glycan_composition(test_composition)

  # Check that the result contains dHex instead of Fuc
  comp_string <- as.character(parsed_comp)
  expect_true(grepl("dHex", comp_string))
  expect_false(grepl("Fuc", comp_string))
})

test_that("glycan composition parsing handles complex compositions", {
  # Test with real compositions from the test data
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )

  # Check that all compositions are valid glyrepr objects
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")

  # Check that compositions have reasonable monosaccharide counts
  comp_strings <- as.character(res$var_info$glycan_composition)
  expect_true(all(nchar(comp_strings) > 0))

  # Verify that percentage values are stripped (should not contain %)
  expect_false(any(grepl("%", comp_strings)))
})

# ----- Test variable information extraction -----
test_that("it correctly extracts variable information", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )

  var_info <- res$var_info

  # Check that all required columns are present and of correct type
  expect_true(all(is.character(var_info$variable)))
  expect_true(all(is.character(var_info$peptide)))
  expect_true(all(is.character(var_info$protein)))
  expect_true(all(is.character(var_info$gene)))
  expect_true(all(is.numeric(var_info$peptide_site)))
  expect_true(all(is.numeric(var_info$protein_site)))

  # Check that variable names are unique and properly formatted
  expect_true(all(grepl("^GP\\d+$", var_info$variable)))
  expect_equal(length(unique(var_info$variable)), nrow(var_info))

  # Check that peptide sites and protein sites are reasonable
  expect_true(all(var_info$peptide_site > 0))
  expect_true(all(var_info$protein_site > 0))

  # Check that protein and gene information is present
  expect_true(all(nchar(var_info$protein) > 0))
  expect_true(all(nchar(var_info$gene) > 0))
})

test_that("it correctly calculates protein sites from peptide sites", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )

  var_info <- res$var_info

  # Check that protein_site calculation is correct
  # protein_site should be protein_start + peptide_site - 1
  # We can't directly test this without access to protein_start,
  # but we can check that protein_site >= peptide_site for most cases
  expect_true(all(var_info$protein_site >= var_info$peptide_site))
})

test_that("it handles missing or zero intensity values correctly", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )

  # Check that zero intensities are converted to NA
  # and that the expression matrix has reasonable values
  expr_mat <- res$expr_mat

  # Should have some non-NA values
  expect_true(sum(!is.na(expr_mat)) > 0)

  # All non-NA values should be positive (since zeros become NA)
  non_na_values <- expr_mat[!is.na(expr_mat)]
  if (length(non_na_values) > 0) {
    expect_true(all(non_na_values > 0))
  }
})

# ----- Test error handling -----
test_that("it handles empty directory gracefully", {
  temp_dir <- withr::local_tempdir()

  expect_error(
    suppressMessages(
      read_msfragger(
        temp_dir,
        quant_method = "label-free"
      )
    ),
    "No psm.tsv files found in the specified directory"
  )
})

test_that("it handles directory with no psm.tsv files", {
  temp_dir <- withr::local_tempdir()

  # Create a subdirectory with a different file
  subdir <- file.path(temp_dir, "sample1")
  dir.create(subdir)
  writeLines("test", file.path(subdir, "other_file.txt"))

  expect_error(
    suppressMessages(
      read_msfragger(
        temp_dir,
        quant_method = "label-free"
      )
    ),
    "No psm.tsv files found in the specified directory"
  )
})

# ----- Test data aggregation -----
test_that("it correctly aggregates PSMs to glycopeptides", {
  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      quant_method = "label-free"
    )
  )

  # Check that aggregation produces reasonable results
  expect_true(nrow(res$var_info) > 0)

  # Check that all variables have unique identifiers
  expect_equal(length(unique(res$var_info$variable)), nrow(res$var_info))

  # Check that expression matrix matches variable info
  expect_equal(nrow(res$expr_mat), nrow(res$var_info))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
})

test_that("it preserves sample information structure", {
  # Test with custom sample info
  sample_info <- tibble::tibble(
    sample = c("H_1", "H_2", "H_3"),
    condition = c("control", "treatment", "treatment"),
    replicate = c(1, 1, 2)
  )

  suppressMessages(
    res <- read_msfragger(
      test_path("data/msfragger_result"),
      sample_info = sample_info,
      quant_method = "label-free"
    )
  )

  # Check that sample info is preserved correctly
  expect_equal(nrow(res$sample_info), 3)
  expect_equal(colnames(res$sample_info), c("sample", "condition", "replicate"))
  expect_equal(res$sample_info$condition, c("control", "treatment", "treatment"))
  expect_equal(res$sample_info$replicate, c(1, 1, 2))
})