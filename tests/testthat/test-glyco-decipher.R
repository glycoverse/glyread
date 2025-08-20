# ----- Basic functionality -----
test_that("it returns correct information (label-free)", {
  suppressMessages(
    res <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free"
    )
  )

  # Check expected columns
  expected_cols <- c("variable", "protein", "protein_site", "glycan_composition")
  expect_true(all(expected_cols %in% colnames(res$var_info)))

  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  expect_type(res$var_info$protein_site, "integer")
  expect_equal(colnames(res$sample_info), c("sample"))
  
  # Check sample names and data structure
  expected_samples <- c("20241224-LXJ-Nglyco-H_1", "20241224-LXJ-Nglyco-H_2", "20241224-LXJ-Nglyco-H_3")
  expect_true(all(expected_samples %in% colnames(res$expr_mat)))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  
  expect_equal(res$meta_data, list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  ))
})


# ----- Data processing -----
test_that("protein inference and site resolution work correctly", {
  suppressMessages(
    res <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free"
    )
  )
  
  # Proteins should be properly inferred (no semicolons, proper format)
  proteins <- res$var_info$protein
  expect_false(any(stringr::str_detect(proteins, ";")))
  expect_true(all(stringr::str_detect(proteins, "^[A-Z0-9]+$")))
  
  # Protein sites should be integers, some may be NA for uncertain sites
  protein_sites <- res$var_info$protein_site
  expect_type(protein_sites, "integer")
  expect_true(any(!is.na(protein_sites)))  # Some sites should be determined
})


test_that("glycan compositions are correctly processed", {
  suppressMessages(
    res <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free"
    )
  )
  
  compositions <- res$var_info$glycan_composition
  expect_s3_class(compositions, "glyrepr_composition")
  
  # Check that Fuc is converted to dHex and adducts are removed
  composition_strings <- as.character(compositions)
  expect_false(any(stringr::str_detect(composition_strings, "Fuc")))
  expect_false(any(stringr::str_detect(composition_strings, "\\+")))
  expect_true(any(stringr::str_detect(composition_strings, "dHex")))
})


test_that("expression matrix and variables are correctly processed", {
  suppressMessages(
    res <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free"
    )
  )

  # Check expression matrix
  expect_true(is.numeric(res$expr_mat))
  expect_equal(nrow(res$expr_mat), nrow(res$var_info))
  expect_equal(ncol(res$expr_mat), nrow(res$sample_info))
  expect_true(all(is.na(res$expr_mat) | res$expr_mat > 0))

  # Check variable naming
  variables <- res$var_info$variable
  expect_true(all(stringr::str_starts(variables, "GP")))
  expect_equal(length(variables), length(unique(variables)))
})


# ----- Sample information handling -----
test_that("it handles sample information correctly", {
  # Default sample info
  suppressMessages(
    res1 <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free"
    )
  )
  expected_samples <- c("20241224-LXJ-Nglyco-H_1", "20241224-LXJ-Nglyco-H_2", "20241224-LXJ-Nglyco-H_3")
  expect_true(all(expected_samples %in% res1$sample_info$sample))
  expect_equal(colnames(res1$sample_info), c("sample"))

  # Custom sample info
  sample_info <- tibble::tibble(
    sample = expected_samples,
    group = c("A", "B", "C")
  )
  suppressMessages(
    res2 <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      sample_info = sample_info,
      quant_method = "label-free"
    )
  )
  expect_equal(colnames(res2$sample_info), c("sample", "group"))
})


test_that("it raises an error for inconsistent sample information", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  expect_error(
    suppressMessages(
      read_glyco_decipher(
        test_path("data/glyco-decipher-result.csv"),
        sample_info = sample_info,
        quant_method = "label-free"
      )
    )
  )
})


test_that("sample name converter works", {
  sample_name_converter <- function(x) {
    stringr::str_replace(x, "20241224-LXJ-Nglyco-H_", "Sample_")
  }

  suppressMessages(
    res <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free",
      sample_name_converter = sample_name_converter
    )
  )

  expected_samples <- c("Sample_1", "Sample_2", "Sample_3")
  expect_true(all(expected_samples %in% colnames(res$expr_mat)))
})


# ----- Parameter validation -----
test_that("it validates parameters correctly", {
  # Invalid file path
  expect_error(
    read_glyco_decipher("nonexistent.csv"),
    "File does not exist"
  )

  # Invalid quant_method
  expect_error(
    read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "invalid_method"
    )
  )

  # TMT not supported
  expect_error(
    read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "TMT"
    ),
    "TMT quantification is not supported yet"
  )
})


test_that("it handles different glycan types and orgdb parameters", {
  # O-linked glycan type
  suppressMessages(
    res1 <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free",
      glycan_type = "O"
    )
  )
  expect_equal(res1$meta_data$glycan_type, "O")

  # Custom orgdb (should work even if package not available)
  suppressMessages(
    res2 <- read_glyco_decipher(
      test_path("data/glyco-decipher-result.csv"),
      quant_method = "label-free",
      orgdb = "org.Mm.eg.db"
    )
  )
  expect_s3_class(res2, "glyexp_experiment")
})