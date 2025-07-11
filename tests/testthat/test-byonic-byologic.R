# ----- Label-free basic functionality -----
test_that("it returns correct information (label-free)", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  # Check expected columns (sample and value are not in var_info after template processing)
  expected_cols <- c(
    "variable", "peptide", "protein", "protein_site",
    "glycan_composition", "peptide_site"
  )
  expect_true(all(expected_cols %in% colnames(res$var_info)))
  
  # Gene column is optional (requires org.Hs.eg.db)
  if ("gene" %in% colnames(res$var_info)) {
    expect_type(res$var_info$gene, "character")
  }

  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  expect_equal(colnames(res$sample_info), c("sample"))
  expect_true(all(c("20241224-LXJ-Nglyco-H_1", "20241224-LXJ-Nglyco-H_2", "20241224-LXJ-Nglyco-H_3") %in% colnames(res$expr_mat)))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(res$meta_data, list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  ))
})


# ----- Multisite glycopeptide filtering -----
test_that("multisite glycopeptides are filtered in full workflow", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )
  
  # Check that no glycan compositions contain commas (indicating multisite)
  compositions <- as.character(res$var_info$glycan_composition)
  expect_false(any(stringr::str_detect(compositions, ",")))
  
  # Verify that the function processes the data and removes multisite entries
  # This is inferred from the fact that all returned compositions are single-site
  expect_true(nrow(res$var_info) > 0)
  expect_true(all(nchar(compositions) > 0))
})


# ----- Row hierarchy collapse -----
test_that("hierarchical rows are correctly collapsed", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )
  
  # Check that we have multiple samples in the expression matrix
  expect_true(ncol(res$expr_mat) > 1)  # Multiple samples

  # Check that each variable appears only once (no duplicates in var_info)
  expect_equal(length(res$var_info$variable), length(unique(res$var_info$variable)))
})


# ----- Protein accession extraction -----
test_that("protein accessions are correctly extracted", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )
  
  # Check protein extraction from UniProt format
  proteins <- res$var_info$protein
  expect_true(all(nchar(proteins) > 0))
  
  # Check some specific examples we know should be in the data (only glycosylated proteins)
  expect_true(any(stringr::str_detect(proteins, "P02765")))  # FETUA_HUMAN (has glycosylation)
  # Note: P02768 (ALBU_HUMAN) is filtered out because it has no glycan modifications in this dataset
})


# ----- Glycan composition processing -----
test_that("glycan compositions are correctly processed and parsed", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )
  
  compositions <- res$var_info$glycan_composition
  expect_s3_class(compositions, "glyrepr_composition")
  
  # Check that Fuc is converted to dHex
  composition_strings <- as.character(compositions)
  expect_false(any(stringr::str_detect(composition_strings, "Fuc")))
  expect_true(any(stringr::str_detect(composition_strings, "dHex")))
  
  # Check that we have expected glycan types
  expect_true(any(stringr::str_detect(composition_strings, "HexNAc")))
  expect_true(any(stringr::str_detect(composition_strings, "Hex")))
  expect_true(any(stringr::str_detect(composition_strings, "NeuAc")))
})


# ----- Peptide processing -----
test_that("peptides are correctly processed", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )
  
  peptides <- res$var_info$peptide
  
  # Check that peptides are uppercase
  expect_true(all(peptides == stringr::str_to_upper(peptides)))
  
  # Check that peptides don't contain flanking amino acids (no dots)
  expect_false(any(stringr::str_detect(peptides, "\\.")))
  
  # Check that we have some expected peptides
  expect_true(any(stringr::str_detect(peptides, "AALAA")))  # From FETUA_HUMAN
})


# ----- Site information -----
test_that("peptide and protein sites are correctly calculated", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )
  
  # Check that peptide_site is integer and reasonable
  expect_type(res$var_info$peptide_site, "integer")
  expect_true(all(res$var_info$peptide_site > 0))
  expect_true(all(res$var_info$peptide_site <= 50))  # Reasonable range for peptide sites
  
  # Check that protein_site is integer and reasonable
  expect_type(res$var_info$protein_site, "integer")
  expect_true(all(res$var_info$protein_site > 0))
  expect_true(all(res$var_info$protein_site <= 1000))  # Reasonable range for protein sites
  
  # Check that protein_site >= peptide_site (since protein_site = start_aa + peptide_site - 1)
  expect_true(all(res$var_info$protein_site >= res$var_info$peptide_site))
})


# ----- Variable naming -----
test_that("variables are correctly named with GP prefix", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  variables <- res$var_info$variable
  expect_true(all(stringr::str_starts(variables, "GP")))
  expect_equal(length(variables), length(unique(variables)))  # All unique
})


# ----- Expression matrix processing -----
test_that("intensity values are correctly processed", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  expr_mat <- res$expr_mat

  # Check that zero values are converted to NA
  expect_true(all(is.na(expr_mat) | expr_mat > 0))

  # Check that the matrix is numeric
  expect_true(is.numeric(expr_mat))

  # Check dimensions
  expect_equal(nrow(expr_mat), nrow(res$var_info))
  expect_equal(ncol(expr_mat), nrow(res$sample_info))
})


# ----- Sample name extraction -----
test_that("sample names are correctly extracted from MS alias names", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  # Check that sample names match the expected pattern
  expected_samples <- c("20241224-LXJ-Nglyco-H_1", "20241224-LXJ-Nglyco-H_2", "20241224-LXJ-Nglyco-H_3")
  expect_true(all(expected_samples %in% res$sample_info$sample))
  expect_true(all(expected_samples %in% colnames(res$expr_mat)))
})


# ----- Error handling -----
test_that("TMT quantification throws error", {
  expect_error(
    read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "TMT"
    ),
    "TMT quantification is not supported yet"
  )
})


test_that("it handles O-linked glycan type", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free",
      glycan_type = "O"
    )
  )
  expect_equal(res$meta_data$glycan_type, "O")
})


# ----- Data consistency -----
test_that("all output data types are consistent", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  # Check that var_info has correct types
  expect_type(res$var_info$variable, "character")
  expect_type(res$var_info$peptide, "character")
  expect_type(res$var_info$protein, "character")
  expect_type(res$var_info$peptide_site, "integer")
  expect_type(res$var_info$protein_site, "integer")
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
})


# ----- Sample name converter -----
test_that("sample name converter works", {
  sample_name_converter <- function(x) {
    stringr::str_replace(x, "20241224-LXJ-Nglyco-H_", "Sample_")
  }

  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free",
      sample_name_converter = sample_name_converter
    )
  )

  expected_samples <- c("Sample_1", "Sample_2", "Sample_3")
  expect_true(all(expected_samples %in% colnames(res$expr_mat)))
  expect_true(all(expected_samples %in% res$sample_info$sample))
})


# ----- Parameter validation -----
test_that("it validates quant_method parameter", {
  expect_error(
    read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "invalid_method"
    ),
    "must be one of"
  )
})


test_that("it validates glycan_type parameter", {
  expect_error(
    read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free",
      glycan_type = "invalid_type"
    ),
    "must be one of"
  )
})


test_that("it validates file path argument", {
  expect_error(
    read_byonic_byologic("nonexistent.csv"),
    "File does not exist"
  )
})


# ----- Sample information handling -----
test_that("it provides a default sample information tibble (label-free)", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  # Should have all three samples from the data
  expected_samples <- c("20241224-LXJ-Nglyco-H_1", "20241224-LXJ-Nglyco-H_2", "20241224-LXJ-Nglyco-H_3")
  expect_true(all(expected_samples %in% res$sample_info$sample))
  expect_equal(colnames(res$sample_info), c("sample"))
})


test_that("it accepts a sample_info tibble", {
  sample_info <- tibble::tibble(
    sample = c("20241224-LXJ-Nglyco-H_1", "20241224-LXJ-Nglyco-H_2", "20241224-LXJ-Nglyco-H_3"),
    group = c("A", "B", "C")
  )

  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      sample_info = sample_info,
      quant_method = "label-free"
    )
  )

  expect_equal(colnames(res$sample_info), c("sample", "group"))
  expect_equal(res$sample_info$group, c("A", "B", "C"))
})


test_that("it raises an error when samples in sample information is inconsistent with result", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  new_sample_info_path <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(sample_info, new_sample_info_path)
  expect_error(
    suppressMessages(
      read_byonic_byologic(
        test_path("data/byonic-byologic-LFQ-result.csv"),
        new_sample_info_path,
        quant_method = "label-free"
      )
    )
  )
})


# ----- Data filtering tests -----
test_that("only glycosylated peptides are included", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  # All entries should have glycan compositions (no NA values)
  expect_false(any(is.na(res$var_info$glycan_composition)))
  expect_true(all(nchar(as.character(res$var_info$glycan_composition)) > 0))
})


# ----- OrgDb parameter -----
test_that("it accepts custom orgdb parameter", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free",
      orgdb = "org.Mm.eg.db"  # Mouse database (won't be installed, but should not error)
    )
  )

  # Should still work even if the orgdb is not available
  expect_s3_class(res, "glyexp_experiment")
})
