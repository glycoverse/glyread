# ----- Label-free -----
test_that("it returns correct information (label-free)", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      test_path("data/pglyco3-LFQ-result-sample-info.csv"),
      quant_method = "label-free"
    )
  )

  expect_equal(
    colnames(res$var_info),
    c(
      "variable", "peptide", "peptide_site", "glycan_composition", "glycan_structure",
      "protein", "gene", "protein_site"
    )
  )
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  expect_s3_class(res$var_info$glycan_structure, "glyrepr_structure")
  expect_type(res$var_info$protein_site, "integer")
  expect_equal(colnames(res$sample_info), c("sample", "group"))
  expect_equal(colnames(res$expr_mat), c("20250315_LiangYuying_GP09_1", "20250315_LiangYuying_GP09_2"))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(res$meta_data, list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  ))
})


test_that("it accepts a sample_info tibble", {
  suppressMessages({
    sample_info <- readr::read_csv(test_path("data/pglyco3-LFQ-result-sample-info.csv"))
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      sample_info = sample_info,
      quant_method = "label-free"
    )
  })
  expect_equal(colnames(res$sample_info), c("sample", "group"))
})


test_that("it accepts a sample_info data.frame", {
  suppressMessages({
    # Create a data.frame (not tibble)
    sample_info <- data.frame(
      sample = c("20250315_LiangYuying_GP09_1", "20250315_LiangYuying_GP09_2"),
      group = c("control", "treatment"),
      stringsAsFactors = FALSE
    )
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      sample_info = sample_info,
      quant_method = "label-free"
    )
  })
  expect_equal(colnames(res$sample_info), c("sample", "group"))
})


test_that("it provides a default sample information tibble (label-free)", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free"
    )
  )
  expect_equal(res$sample_info, tibble::tibble(sample = c("20250315_LiangYuying_GP09_2", "20250315_LiangYuying_GP09_1")))
})


test_that("it raises an error when samples in sample information is inconsistent with pGlyco3 result (label-free)", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  new_sample_info_path <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(sample_info, new_sample_info_path)
  expect_error(
    suppressMessages(
      read_pglyco3(
        test_path("data/pglyco3-LFQ-result.txt"),
        new_sample_info_path,
        quant_method = "label-free"
      )
    )
  )
})


test_that("it works with O-glycan type", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free",
      glycan_type = "O"
    )
  )
  expect_equal(res$meta_data$glycan_type, "O")
})


test_that("it works with sample_name_converter", {
  converter <- function(x) {
    stringr::str_replace(x, "20250315_LiangYuying_GP09_", "Sample_")
  }

  sample_info <- tibble::tibble(
    sample = c("Sample_2", "Sample_1"),  # Note: order matches the actual sample order
    group = c("A", "B")
  )

  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      sample_info = sample_info,
      quant_method = "label-free",
      sample_name_converter = converter
    )
  )

  expect_equal(colnames(res$expr_mat), c("Sample_2", "Sample_1"))
  expect_equal(res$sample_info$sample, c("Sample_2", "Sample_1"))
})


# ----- TMT -----
test_that("it raises an error for TMT quantification (not supported yet)", {
  expect_error(
    suppressMessages(
      read_pglyco3(
        test_path("data/pglyco3-LFQ-result.txt"),
        quant_method = "TMT"
      )
    ),
    "TMT quantification is not supported yet"
  )
})


# ----- Argument validation -----
test_that("it validates file path argument", {
  expect_error(
    read_pglyco3("nonexistent.txt"),
    "File does not exist"
  )
})


test_that("it validates quant_method argument", {
  expect_error(
    read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "invalid"
    )
  )
})


test_that("it validates glycan_type argument", {
  expect_error(
    read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      glycan_type = "invalid"
    )
  )
})


test_that("it validates sample_name_converter argument", {
  expect_error(
    read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      sample_name_converter = "not_a_function"
    )
  )
})


# ----- Data structure tests -----
test_that("glycan composition parsing works correctly", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free"
    )
  )
  
  # Check that all glycan compositions are parsed
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  compositions <- res$var_info$glycan_composition
  expect_true(length(compositions) > 0)
})


test_that("glycan structure parsing works correctly", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free"
    )
  )
  
  # Check that all glycan structures are parsed
  expect_s3_class(res$var_info$glycan_structure, "glyrepr_structure")
  structures <- res$var_info$glycan_structure
  expect_true(length(structures) > 0)
})


test_that("protein and gene information is correctly processed", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free"
    )
  )
  
  genes <- res$var_info$gene
  # Check that gene names are properly processed (after protein inference)
  expect_true(all(nchar(genes) > 0))
  expect_true(all(!is.na(genes)))
})


test_that("expression matrix has correct dimensions and values", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free"
    )
  )
  
  # Check dimensions
  expect_equal(nrow(res$expr_mat), nrow(res$var_info))
  expect_equal(ncol(res$expr_mat), nrow(res$sample_info))
  
  # Check that values are numeric and non-negative
  expect_type(res$expr_mat, "double")
  expect_true(all(res$expr_mat >= 0, na.rm = TRUE))
})


test_that("variable identifiers are unique", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free"
    )
  )
  
  variables <- res$var_info$variable
  expect_equal(length(variables), length(unique(variables)))
  expect_true(all(stringr::str_starts(variables, "GP")))
})

# ----- Protein inference tests -----
test_that("it performs protein inference with parsimony method by default", {
  suppressMessages(
    res <- read_pglyco3(
      test_path("data/pglyco3-LFQ-result.txt"),
      quant_method = "label-free"
    )
  )
  
  # Check that protein inference columns exist
  expect_true("protein" %in% colnames(res$var_info))
  expect_true("gene" %in% colnames(res$var_info))
  expect_true("protein_site" %in% colnames(res$var_info))
  
  # Check that original columns are removed
  expect_false("proteins" %in% colnames(res$var_info))
  expect_false("genes" %in% colnames(res$var_info))
  expect_false("protein_sites" %in% colnames(res$var_info))
  
  # Check that protein_site is integer
  expect_type(res$var_info$protein_site, "integer")
})


# ----- PSM aggregation tests -----
test_that("PSMs are correctly aggregated to glycopeptides", {
  # Read data and double the rows to simulate multiple PSMs per glycopeptide
  suppressMessages(df <- readr::read_tsv(test_path("data/pglyco3-LFQ-result.txt")))
  new_df <- dplyr::bind_rows(df, df)
  withr::with_dir(tempdir(), {
    readr::write_tsv(new_df, "test.txt")
    suppressMessages(
      res <- read_pglyco3("test.txt")
    )
  })

  # Check that aggregation worked correctly - should have same number of unique glycopeptides
  # The original data has 6 unique glycopeptides after aggregation
  expect_equal(nrow(res$var_info), 6)
})
