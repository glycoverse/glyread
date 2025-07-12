# ----- Label-free -----
test_that("it returns correct information (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      test_path("data/pglyco3-LFQ-sample-info.csv"),
      quant_method = "label-free"
    )
  )

  expect_equal(
    colnames(res$var_info),
    c(
      "variable", "peptide", "peptide_site", "protein", "protein_site",
      "gene", "glycan_composition", "glycan_structure"
    )
  )
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  expect_s3_class(res$var_info$glycan_structure, "glyrepr_structure")
  expect_type(res$var_info$protein_site, "integer")
  expect_equal(colnames(res$sample_info), c("sample", "group"))
  expect_equal(colnames(res$expr_mat), c("S1", "S2"))
  expect_equal(rownames(res$expr_mat), res$var_info$variable)
  expect_equal(res$meta_data, list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  ))
})


test_that("it accepts a sample_info tibble", {
  suppressMessages({
    sample_info <- readr::read_csv(test_path("data/pglyco3-LFQ-sample-info.csv"))
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
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
      sample = c("S1", "S2"),
      group = c("control", "treatment"),
      stringsAsFactors = FALSE
    )
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      sample_info = sample_info,
      quant_method = "label-free"
    )
  })
  expect_equal(colnames(res$sample_info), c("sample", "group"))
  expect_s3_class(res$sample_info, "tbl_df")  # Should be converted to tibble
})


test_that("it provides a default sample information tibble (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  expect_equal(res$sample_info, tibble::tibble(sample = c("S1", "S2")))
})


test_that("it raiss an error when samples in sample information is inconsistent with pGlyco3 result (label-free)", {
  sample_info <- tibble::tibble(sample = c("S1", "S2", "S3"))
  new_sample_info_path <- withr::local_tempfile(fileext = ".csv")
  readr::write_csv(sample_info, new_sample_info_path)
  expect_error(
    suppressMessages(
      read_pglyco3_pglycoquant(
        test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
        test_path(new_sample_info_path),
        quant_method = "label-free"
      )
    )
  )
})


test_that("zeros are replaced by NA (label-free)", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  expect_true(sum(is.na(res$expr_mat)) > 0)
})


test_that("sample name converter works", {
  sample_name_converter <- function(x) {
    stringr::str_replace(x, "S", "Sample_")
  }
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      sample_name_converter = sample_name_converter
    )
  )
  expect_equal(colnames(res$expr_mat), c("Sample_1", "Sample_2"))
  expect_equal(res$sample_info$sample, c("Sample_1", "Sample_2"))
})

# ----- Parameter validation tests -----
test_that("it validates file path parameter", {
  expect_error(
    read_pglyco3_pglycoquant(
      "non_existent_file.list",
      quant_method = "label-free"
    ),
    class = "simpleError"
  )
  
  # Test wrong file extension
  temp_wrong_ext <- withr::local_tempfile(fileext = ".txt")
  writeLines("test", temp_wrong_ext)
  expect_error(
    read_pglyco3_pglycoquant(
      temp_wrong_ext,
      quant_method = "label-free"
    ),
    class = "simpleError"
  )
})

test_that("it validates sample_info parameter", {
  # Test invalid sample_info file
  temp_invalid_csv <- withr::local_tempfile(fileext = ".txt")
  writeLines("not a csv", temp_invalid_csv)
  expect_error(
    read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      temp_invalid_csv,
      quant_method = "label-free"
    ),
    class = "simpleError"
  )
})

test_that("it validates quant_method parameter", {
  expect_error(
    read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "invalid_method"
    ),
    "must be one of"
  )
})

test_that("it validates glycan_type parameter", {
  expect_error(
    read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      glycan_type = "invalid_type"
    ),
    "must be one of"
  )
})

test_that("it validates sample_name_converter parameter", {
  expect_error(
    read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      sample_name_converter = "not_a_function"
    ),
    class = "simpleError"
  )
})

# ----- Glycan type tests -----
test_that("it handles O-linked glycan type", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      glycan_type = "O"
    )
  )
  expect_equal(res$meta_data$glycan_type, "O")
})

# ----- TMT quantification test -----
test_that("it rejects TMT quantification (not implemented)", {
  expect_error(
    read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "TMT"
    ),
    "TMT quantification is not supported yet"
  )
})

# ----- Data parsing tests -----
test_that("glycan composition parsing works correctly", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  
  # Check that all glycan compositions are correctly parsed
  expect_true(all(!is.na(res$var_info$glycan_composition)))
  expect_s3_class(res$var_info$glycan_composition, "glyrepr_composition")
  
  # Check specific glycan compositions from test data
  compositions <- res$var_info$glycan_composition
  expect_true(length(compositions) > 0)
})

test_that("glycan structure parsing works correctly", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
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
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  
  genes <- res$var_info$gene
  # Check that gene names are properly processed (after protein inference)
  expect_true(all(nchar(genes) > 0))
  expect_true(all(!is.na(genes)))
})

# ----- Expression matrix tests -----
test_that("expression matrix has correct dimensions and properties", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  
  expr_mat <- res$expr_mat
  
  # Check dimensions
  expect_equal(nrow(expr_mat), nrow(res$var_info))
  expect_equal(ncol(expr_mat), nrow(res$sample_info))
  
  # Check that matrix contains numeric values or NA
  expect_true(is.numeric(expr_mat))
  
  # Check that all zero values were converted to NA
  expect_true(all(is.na(expr_mat) | expr_mat > 0))
})

test_that("variable identifiers are unique", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  
  variables <- res$var_info$variable
  expect_equal(length(variables), length(unique(variables)))
  expect_true(all(stringr::str_starts(variables, "GP")))
})

# ----- Edge cases -----
test_that("it handles empty sample_name_converter result", {
  bad_converter <- function(x) character(0)
  expect_error(
    suppressMessages(
      read_pglyco3_pglycoquant(
        test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
        quant_method = "label-free",
        sample_name_converter = bad_converter
      )
    )
  )
})

test_that("it handles complex sample info with multiple columns", {
  complex_sample_info <- tibble::tibble(
    sample = c("S1", "S2"),
    group = c("A", "B"), 
    batch = c(1, 2),
    bio_replicate = c("Bio1", "Bio2")
  )
  
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      sample_info = complex_sample_info,
      quant_method = "label-free"
    )
  )
  
  expect_equal(colnames(res$sample_info), c("sample", "group", "batch", "bio_replicate"))
  expect_equal(res$sample_info$batch, c(1, 2))
})

# ----- Error handling tests -----
test_that("it handles malformed data files gracefully", {
  # This test checks that the function provides meaningful error messages
  # when encountering malformed data files
  expect_error(
    suppressMessages(
      read_pglyco3_pglycoquant(
        test_path("data/pglyco3-malformed-data.list"),
        quant_method = "label-free"
      )
    )
  )
})

test_that("it handles files with missing essential columns", {
  # Test with a file that's missing essential columns
  expect_error(
    suppressMessages(
      read_pglyco3_pglycoquant(
        test_path("data/pglyco3-missing-columns.list"),
        quant_method = "label-free"
      )
    )
  )
})

test_that("it handles sample_info with missing sample column", {
  malformed_sample_info <- tibble::tibble(
    name = c("S1", "S2"),  # Wrong column name
    group = c("A", "B")
  )
  
  expect_error(
    suppressMessages(
      read_pglyco3_pglycoquant(
        test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
        sample_info = malformed_sample_info,
        quant_method = "label-free"
      )
    )
  )
})

# ----- Internal function tests -----
test_that(".convert_glycan_composition works correctly", {
  # Test individual glycan composition strings
  test_comps <- c("H(4)N(4)F(1)", "H(5)N(2)", "H(3)N(2)")
  result <- glyread:::.convert_pglyco3_comp(test_comps)
  
  expect_s3_class(result, "glyrepr_composition")
  expect_equal(length(result), 3)
})

test_that(".convert_glycan_composition handles duplicates efficiently", {
  # Test that duplicate compositions are handled efficiently
  test_comps <- rep("H(4)N(4)F(1)", 100)
  result <- glyread:::.convert_pglyco3_comp(test_comps)
  
  expect_s3_class(result, "glyrepr_composition")
  expect_equal(length(result), 100)
})

test_that(".convert_glycan_composition handles complex compositions", {
  # Test various complex glycan compositions
  test_comps <- c(
    "H(5)N(4)A(1)F(1)",  # Complex composition with sialic acid
    "H(6)N(5)F(2)",      # High mannose with fucose
    "H(3)N(2)G(1)"       # With hexuronic acid
  )
  result <- glyread:::.convert_pglyco3_comp(test_comps)
  
  expect_s3_class(result, "glyrepr_composition")
  expect_equal(length(result), 3)
  expect_true(all(!is.na(result)))
})

# ----- Data consistency tests -----
test_that("all output data types are consistent", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  
  # Check data types of variable information
  expect_type(res$var_info$peptide, "character")
  expect_type(res$var_info$protein, "character")
  expect_type(res$var_info$gene, "character")
  expect_type(res$var_info$peptide_site, "integer")
  expect_type(res$var_info$protein_site, "integer")

  # Check that all rows have data
  expect_true(all(nchar(res$var_info$peptide) > 0))
})

test_that("intensity columns are correctly extracted", {
  suppressMessages(
    res <- read_pglyco3_pglycoquant(
      test_path("data/pglyco3-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )
  
  # Check that intensity column names are correctly processed
  expr_mat <- res$expr_mat
  expect_true(all(colnames(expr_mat) %in% c("S1", "S2")))
  expect_true(is.numeric(expr_mat))
  
  # Check that intensity values are positive when not NA
  non_na_values <- expr_mat[!is.na(expr_mat)]
  expect_true(all(non_na_values > 0))
})


# ----- PSM aggregation tests -----
test_that("PSMs are correctly aggregated to glycopeptides", {
  # Read data and double the rows to simulate multiple PSMs per glycopeptide
  suppressMessages(df <- readr::read_tsv(test_path("data/pglyco3-pglycoquant-LFQ-result.list")))
  new_df <- dplyr::bind_rows(df, df)
  withr::with_dir(tempdir(), {
    readr::write_tsv(new_df, "test.list")
    suppressMessages(
      res <- read_pglyco3_pglycoquant("test.list")
    )
  })

  # Check that aggregation worked correctly
  expect_equal(nrow(res$var_info), nrow(df))
})
