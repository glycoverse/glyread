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
    "variable",
    "peptide",
    "protein",
    "protein_site",
    "glycan_composition",
    "peptide_site"
  )
  expect_true(all(expected_cols %in% colnames(res$var_info)))

  # Gene column is optional (requires org.Hs.eg.db)
  if ("gene" %in% colnames(res$var_info)) {
    expect_type(res$var_info$gene, "character")
  }

  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  expect_equal(colnames(res$sample_info), c("sample"))
  expect_true(all(
    c(
      "20241224-LXJ-Nglyco-H_1",
      "20241224-LXJ-Nglyco-H_2",
      "20241224-LXJ-Nglyco-H_3"
    ) %in%
      colnames(res$expr_mat)
  ))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(
    res$meta_data,
    list(
      exp_type = "glycoproteomics",
      glycan_type = "N",
      quant_method = "label-free"
    )
  )
})


# ----- Multisite glycopeptide handling -----
test_that("multisite glycopeptides are handled correctly in full workflow", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  # All entries should be present (multisite glycopeptides are no longer filtered out)
  expect_true(nrow(res$var_info) > 0)

  # Multisite glycopeptides should be expanded before composition parsing,
  # so parsed output should not retain comma-separated compositions.
  compositions <- as.character(res$var_info$glycan_composition)
  multisite_entries <- stringr::str_detect(compositions, ",")
  expect_false(any(multisite_entries))

  # Expanded entries should have valid protein_site values.
  expect_true(all(!is.na(res$var_info$protein_site)))
})


test_that("multisite glycopeptides are expanded into site-specific rows", {
  byologic_path <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(
    tibble::tibble(
      row_number = "1",
      protein_name = "sp|P02765|FETUA_HUMAN",
      sequence = "K.ABCDnEFGHnIK.R",
      glycans = "HexNAc(1)Fuc(1),HexNAc(4)Hex(5)Fuc(1)NeuAc(1)",
      xic_area_summed = 100,
      ms_alias_name = "Sample1",
      mod_summary = "N5(NGlycan/349.1373); N9(NGlycan/2059.7349)",
      start_aa = 10L
    ),
    byologic_path
  )

  suppressMessages(
    res <- read_byonic_byologic(
      byologic_path,
      quant_method = "label-free",
      orgdb = "missing.OrgDb"
    )
  )

  res$var_info <- dplyr::arrange(res$var_info, .data$protein_site)

  expect_equal(nrow(res$var_info), 2L)
  expect_true("gp_id" %in% colnames(res$var_info))
  expect_equal(length(unique(res$var_info$gp_id)), 1L)
  expect_equal(res$var_info$peptide_site, c(5L, 9L))
  expect_equal(res$var_info$protein_site, c(14L, 18L))
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  expect_equal(
    as.character(res$var_info$glycan_composition),
    c("HexNAc(1)dHex(1)", "Hex(5)HexNAc(4)dHex(1)NeuAc(1)")
  )
  expect_equal(as.numeric(res$expr_mat[, "Sample1"]), c(100, 100))
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
  expect_true(ncol(res$expr_mat) > 1) # Multiple samples

  # Check that each variable appears only once (no duplicates in var_info)
  expect_equal(
    length(res$var_info$variable),
    length(unique(res$var_info$variable))
  )
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
  expect_true(any(stringr::str_detect(proteins, "P02765"))) # FETUA_HUMAN (has glycosylation)
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
  expect_true(any(stringr::str_detect(peptides, "AALAA"))) # From FETUA_HUMAN
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
  expect_true(all(res$var_info$peptide_site <= 50)) # Reasonable range for peptide sites

  # Check that protein_site is integer and reasonable
  expect_type(res$var_info$protein_site, "integer")
  expect_true(all(res$var_info$protein_site > 0))
  expect_true(all(res$var_info$protein_site <= 1000)) # Reasonable range for protein sites

  # Check that protein_site >= peptide_site (since protein_site = start_aa + peptide_site - 1)
  expect_true(all(res$var_info$protein_site >= res$var_info$peptide_site))
})


# ----- Variable naming -----
test_that("variables are correctly named with meaningful IDs", {
  suppressMessages(
    res <- read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free"
    )
  )

  variables <- res$var_info$variable
  # Variables should be meaningful: protein-site-glycan pattern
  expect_true(all(stringr::str_detect(variables, ".+?-\\d+-.+")))
  # And contain glycan composition info
  expect_true(all(stringr::str_detect(variables, "Hex")))
  expect_true(all(stringr::str_detect(variables, "HexNAc")))
  expect_equal(length(variables), length(unique(variables))) # All unique
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
  expected_samples <- c(
    "20241224-LXJ-Nglyco-H_1",
    "20241224-LXJ-Nglyco-H_2",
    "20241224-LXJ-Nglyco-H_3"
  )
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
      glycan_type = "O-GalNAc"
    )
  )
  expect_equal(res$meta_data$glycan_type, "O-GalNAc")
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
    " Must be element of set"
  )
})


test_that("it validates glycan_type parameter", {
  expect_error(
    read_byonic_byologic(
      test_path("data/byonic-byologic-LFQ-result.csv"),
      quant_method = "label-free",
      glycan_type = "invalid_type"
    ),
    " Must be element of set"
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
  expected_samples <- c(
    "20241224-LXJ-Nglyco-H_1",
    "20241224-LXJ-Nglyco-H_2",
    "20241224-LXJ-Nglyco-H_3"
  )
  expect_true(all(expected_samples %in% res$sample_info$sample))
  expect_equal(colnames(res$sample_info), c("sample"))
})


test_that("it accepts a sample_info tibble", {
  sample_info <- tibble::tibble(
    sample = c(
      "20241224-LXJ-Nglyco-H_1",
      "20241224-LXJ-Nglyco-H_2",
      "20241224-LXJ-Nglyco-H_3"
    ),
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
  expect_equal(res$sample_info$group, factor(c("A", "B", "C")))
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
      orgdb = "org.Mm.eg.db" # Mouse database (won't be installed, but should not error)
    )
  )

  # Should still work even if the orgdb is not available
  expect_s3_class(res, "glyexp_experiment")
})
