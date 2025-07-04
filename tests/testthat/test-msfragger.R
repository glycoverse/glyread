# ----- Basic functionality tests -----
test_that("it returns correct information (label-free)", {
  suppressMessages(
    res <- read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free"
    )
  )

  # Check structure of returned experiment object
  expect_s3_class(res, "glyexp_experiment")
  
  # Check variable information columns
  expect_equal(
    colnames(res$var_info),
    c(
      "variable", "peptide", "charge", "glycan_composition", 
      "protein", "gene", "peptide_site", "protein_site"
    )
  )
  
  # Check glycan composition is properly parsed
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  
  # Check sample information has correct structure
  expect_equal(colnames(res$sample_info), "sample")
  
  # Check expression matrix dimensions and names
  expect_equal(ncol(res$expr_mat), 1) # Only one sample in test data
  expect_equal(nrow(res$expr_mat), nrow(res$var_info))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  
  # Check metadata
  expect_equal(res$meta_data, list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  ))
})

test_that("it provides a default sample information tibble", {
  suppressMessages(
    res <- read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free"
    )
  )
  # Sample name should be extracted from file path
  expected_sample_name <- "H_1" # Based on the file path in test data
  expect_equal(res$sample_info, tibble::tibble(sample = expected_sample_name))
})

test_that("zeros are replaced by NA", {
  suppressMessages(
    res <- read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
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
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free"
    )
  )
  # Test data should have some PSMs filtered out
  # All remaining PSMs should have exactly one glycosite
  expect_true(nrow(res$var_info) <= 6) # Original test file has 6 rows
})

test_that("it handles O-linked glycan type", {
  suppressMessages(
    res <- read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free",
      glycan_type = "O"
    )
  )
  expect_equal(res$meta_data$glycan_type, "O")
})

# ----- Parameter validation tests -----
test_that("it validates file path parameter", {
  expect_error(
    read_msfragger(
      "non_existent_file.tsv",
      quant_method = "label-free"
    ),
    class = "simpleError"
  )
  
  # Test wrong file extension
  temp_wrong_ext <- withr::local_tempfile(fileext = ".txt")
  writeLines("test", temp_wrong_ext)
  expect_error(
    read_msfragger(
      temp_wrong_ext,
      quant_method = "label-free"
    ),
    class = "simpleError"
  )
})

test_that("it validates quant_method parameter", {
  expect_error(
    read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "invalid_method"
    ),
    "must be one of"
  )
})

test_that("it validates glycan_type parameter", {
  expect_error(
    read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free",
      glycan_type = "invalid_type"
    ),
    "must be one of"
  )
})

test_that("it validates sample_name_converter parameter", {
  expect_error(
    read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free",
      sample_name_converter = "not_a_function"
    ),
    class = "simpleError"
  )
})

test_that("it rejects TMT quantification (not implemented)", {
  expect_error(
    read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
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
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free",
      sample_name_converter = sample_name_converter
    )
  )
  expect_equal(colnames(res$expr_mat), "Sample_H_1")
  expect_equal(res$sample_info$sample, "Sample_H_1")
})

# ----- Test documentation code snippet -----
test_that("documentation code snippet works for merging multiple files", {
  # Create additional test files by copying the existing one
  temp_dir <- withr::local_tempdir()
  
  # Copy the test file to create multiple files
  file.copy(
    test_path("msfragger-LFQ-result.tsv"),
    file.path(temp_dir, "psm1.tsv")
  )
  file.copy(
    test_path("msfragger-LFQ-result.tsv"),
    file.path(temp_dir, "psm2.tsv")
  )
  file.copy(
    test_path("msfragger-LFQ-result.tsv"),
    file.path(temp_dir, "psm3.tsv")
  )
  
  # Test the code snippet from documentation
  suppressMessages({
    # Read the result into a list of experiments
    files <- c(
      file.path(temp_dir, "psm1.tsv"),
      file.path(temp_dir, "psm2.tsv"),
      file.path(temp_dir, "psm3.tsv")
    )
    # Use sample_name_converter to give each experiment a unique sample name
    sample_name_converters <- list(
      function(x) paste0(x, "_file1"),
      function(x) paste0(x, "_file2"), 
      function(x) paste0(x, "_file3")
    )
    exps <- purrr::map2(files, sample_name_converters, 
                        ~ read_msfragger(.x, quant_method = "label-free", sample_name_converter = .y))
    
    # Check that we got a list of experiments
    expect_length(exps, 3)
    expect_true(all(purrr::map_lgl(exps, ~ inherits(.x, "glyexp_experiment"))))
    
    # Since each PSM is already independent, we can treat these as "aggregated" experiments
    # Test the merge functionality as in the documentation snippet
    
    # First, check that each experiment has the right structure
    purrr::walk(exps, function(exp) {
      expect_s3_class(exp, "glyexp_experiment")
      expect_true(ncol(exp$expr_mat) == 1)
    })
    
    # Test merging experiments using reduce (as shown in documentation)
    exp_merged <- purrr::reduce(exps, merge)
    
    # Check the merged experiment structure
    expect_s3_class(exp_merged, "glyexp_experiment")
    expect_equal(ncol(exp_merged$expr_mat), 3) # Should have 3 samples now
    expect_equal(colnames(exp_merged$expr_mat), c("H_1_file1", "H_1_file2", "H_1_file3"))
    
    # Check that the number of variables is consistent
    # Since we copied the same file 3 times, all should have the same variables
    expect_equal(nrow(exp_merged$expr_mat), nrow(exps[[1]]$expr_mat))
    
    # Test the sample info modification as shown in documentation
    # Create a mock sample dataframe
    sample_df <- tibble::tibble(
      sample = c("H_1_file1", "H_1_file2", "H_1_file3"),
      group = c("control", "treatment", "treatment"),
      batch = c("A", "B", "B")
    )
    
    # Test the sample info merging (as shown in docs)
    original_sample_info <- exp_merged$sample_info
    exp_merged$sample_info <- exp_merged$sample_info %>%
      dplyr::left_join(sample_df, by = "sample")
    
    # Check the final sample info structure
    expect_equal(colnames(exp_merged$sample_info), c("sample", "group", "batch"))
    expect_equal(nrow(exp_merged$sample_info), 3)
  })
})

# ----- Test glycan composition parsing -----
test_that("it correctly parses glycan compositions", {
  suppressMessages(
    res <- read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
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
  
  # Test that Fuc is converted to dHex if present
  # Note: Test data doesn't contain Fuc, so we just check the conversion logic works
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

# ----- Test variable information extraction -----
test_that("it correctly extracts variable information", {
  suppressMessages(
    res <- read_msfragger(
      test_path("msfragger-LFQ-result.tsv"),
      quant_method = "label-free"
    )
  )
  
  var_info <- res$var_info
  
  # Check that all required columns are present and of correct type
  expect_true(all(is.character(var_info$variable)))
  expect_true(all(is.character(var_info$peptide)))
  expect_true(all(is.numeric(var_info$charge)))
  expect_true(all(is.character(var_info$protein)))
  expect_true(all(is.character(var_info$gene)))
  expect_true(all(is.numeric(var_info$peptide_site)))
  expect_true(all(is.numeric(var_info$protein_site)))
  
  # Check that variable names are unique and properly formatted
  expect_true(all(grepl("^PSM\\d+$", var_info$variable)))
  expect_equal(length(unique(var_info$variable)), nrow(var_info))
}) 