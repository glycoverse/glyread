# ----- Label-free basic functionality -----
test_that("it returns correct information (label-free)", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  # Check expected columns
  expected_cols <- c(
    "variable",
    "peptide",
    "protein",
    "protein_site",
    "glycan_composition",
    "peptide_site"
  )
  expect_true(all(expected_cols %in% colnames(.test_var_info(res))))

  # Gene column is optional (requires clusterProfiler and org.Hs.eg.db)
  if ("gene" %in% colnames(.test_var_info(res))) {
    expect_type(.test_var_info(res)$gene, "character")
  }

  expect_s3_class(.test_var_info(res)$glycan_composition, "glyrepr_composition")
  expect_equal(colnames(.test_sample_info(res)), c("sample"))
  expect_equal(colnames(.test_expr_mat(res)), c("20241224-lxj-nglyco-h_1"))
  expect_equal(rownames(.test_expr_mat(res)), .test_var_info(res)$variable)
  expect_equal(
    .test_metadata(res),
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
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  # All entries should be present (multisite glycopeptides are no longer filtered out)
  expect_true(nrow(.test_var_info(res)) > 0)

  # Multisite glycopeptides should be expanded before composition parsing,
  # so parsed output should not retain comma-separated compositions.
  compositions <- as.character(.test_var_info(res)$glycan_composition)
  multisite_entries <- stringr::str_detect(compositions, ",")
  expect_false(any(multisite_entries))

  # Expanded entries should have valid protein_site values.
  expect_true(all(!is.na(.test_var_info(res)$protein_site)))
})

test_that("multisite Byonic-pGlycoQuant glycopeptides are expanded into site-specific rows", {
  pglycoquant_path <- withr::local_tempfile(fileext = ".list")
  readr::write_tsv(
    tibble::tibble(
      Peptide = "K.NLFLN[+349.13728]HSEN[+2059.73493]ATAK.D",
      `Protein Name` = ">tr|H0Y300|H0Y300_HUMAN Haptoglobin OS=Homo sapiens OX=9606 GN=HP PE=1 SV=4",
      Position = 239L,
      Composition = "HexNAc(1)Fuc(1),HexNAc(4)Hex(5)Fuc(1)NeuAc(1)",
      `Intensity(Sample1)` = 100
    ),
    pglycoquant_path
  )

  suppressMessages(
    res <- read_byonic_pglycoquant(
      pglycoquant_path,
      quant_method = "label-free",
      orgdb = "missing.OrgDb"
    )
  )

  var_info <- dplyr::arrange(.test_var_info(res), .data$protein_site)

  expect_equal(nrow(var_info), 2L)
  expect_false("gp_id" %in% colnames(var_info))
  expect_equal(var_info$peptide_site, c(5L, 9L))
  expect_equal(var_info$protein_site, c(239L, 243L))
  expect_s3_class(var_info$glycan_composition, "glyrepr_composition")
  expect_equal(
    as.character(var_info$glycan_composition),
    c("HexNAc(1)dHex(1)", "Hex(5)HexNAc(4)dHex(1)NeuAc(1)")
  )
  expect_equal(as.numeric(.test_expr_mat(res)[, "Sample1"]), c(100, 100))
})

test_that("multisite Byonic-pGlycoQuant glycopeptides can be dropped", {
  pglycoquant_path <- withr::local_tempfile(fileext = ".list")
  readr::write_tsv(
    tibble::tibble(
      Peptide = c(
        "K.NLFLN[+349.13728]HSEN[+2059.73493]ATAK.D",
        "K.N[+203.07937]GTR.G"
      ),
      `Protein Name` = c(
        ">tr|H0Y300|H0Y300_HUMAN Haptoglobin OS=Homo sapiens OX=9606 GN=HP PE=1 SV=4",
        ">tr|A0A024R6P0|A0A024R6P0_HUMAN Serpin peptidase inhibitor, clade A, member 3 OS=Homo sapiens OX=9606 GN=SERPINA3 PE=3 SV=1"
      ),
      Position = c(239L, 186L),
      Composition = c(
        "HexNAc(1)Fuc(1),HexNAc(4)Hex(5)Fuc(1)NeuAc(1)",
        "HexNAc(1)"
      ),
      `Intensity(Sample1)` = c(100, 50)
    ),
    pglycoquant_path
  )

  suppressMessages(
    res <- read_byonic_pglycoquant(
      pglycoquant_path,
      quant_method = "label-free",
      orgdb = "missing.OrgDb",
      multisite = "drop"
    )
  )

  expect_equal(nrow(.test_var_info(res)), 1L)
  expect_equal(.test_var_info(res)$peptide_site, 1L)
  expect_equal(.test_var_info(res)$protein_site, 186L)
  expect_equal(
    as.character(.test_var_info(res)$glycan_composition),
    "HexNAc(1)"
  )
  expect_equal(as.numeric(.test_expr_mat(res)[, "Sample1"]), 50)
})

test_that("Byonic-pGlycoQuant rows without glycosite markers fail clearly", {
  malformed_df <- tibble::tibble(
    Peptide = "K.NLFLNHSENATAK.D",
    `Protein Name` = ">tr|H0Y300|H0Y300_HUMAN Haptoglobin OS=Homo sapiens OX=9606 GN=HP PE=1 SV=4",
    Position = 239L,
    Composition = "HexNAc(1)Fuc(1),HexNAc(4)Hex(5)Fuc(1)NeuAc(1)",
    `Intensity(Sample1)` = 100
  )

  expect_error(
    glyread:::.expand_byonic_pglycoquant_multisite_rows(malformed_df),
    "Cannot pair Byonic-pGlycoQuant glycan compositions with glycosylation sites"
  )
})

# ----- Peptide site calculation -----
test_that("peptide_site and protein_site are correctly calculated", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  # Check that peptide_site column exists and contains valid integers
  expect_true("peptide_site" %in% colnames(.test_var_info(res)))
  expect_type(.test_var_info(res)$peptide_site, "integer")
  expect_true(all(.test_var_info(res)$peptide_site > 0))

  # Check that protein_site column exists and contains integers
  expect_true("protein_site" %in% colnames(.test_var_info(res)))
  expect_type(.test_var_info(res)$protein_site, "integer")
  expect_true(all(.test_var_info(res)$protein_site > 0))

  # Check that peptide sequences are properly processed (no modification marks)
  peptides <- .test_var_info(res)$peptide
  expect_true(all(nchar(peptides) > 0))
  expect_false(any(stringr::str_detect(peptides, "\\["))) # No modification brackets
  expect_false(any(stringr::str_detect(peptides, "\\."))) # No prefix/suffix dots
})

test_that("protein accessions are correctly extracted", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  # Check protein extraction from UniProt format
  proteins <- .test_var_info(res)$protein
  expect_true(all(nchar(proteins) > 0))

  # Check some specific examples we know should be in the data
  expect_true(any(stringr::str_detect(proteins, "A0A024R6P0")))
  expect_true(any(stringr::str_detect(proteins, "P19652")))
})

# ----- Glycan composition processing -----
test_that("glycan compositions are correctly processed and parsed", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  # Check that all glycan compositions are valid glyrepr objects
  expect_s3_class(.test_var_info(res)$glycan_composition, "glyrepr_composition")
  expect_true(all(!is.na(.test_var_info(res)$glycan_composition)))

  # Check that Fuc has been converted to dHex in the compositions
  # This is verified by checking that the compositions are properly parsed by glyrepr
  # which expects dHex, not Fuc notation
  compositions_str <- as.character(.test_var_info(res)$glycan_composition)
  expect_true(all(nchar(compositions_str) > 0))

  # Verify that standard monosaccharide nomenclature is used
  expect_false(any(stringr::str_detect(compositions_str, "Fuc"))) # Should be dHex
  expect_true(any(stringr::str_detect(compositions_str, "dHex|HexNAc|Hex"))) # Contains standard notation
})

# ----- Variable naming -----
test_that("variables are correctly named with meaningful IDs", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  variables <- .test_var_info(res)$variable
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
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  expr_mat <- .test_expr_mat(res)

  # Check that zero values are converted to NA
  expect_true(all(is.na(expr_mat) | expr_mat > 0))

  # Check that the matrix is numeric
  expect_true(is.numeric(expr_mat))

  # Check dimensions
  expect_equal(nrow(expr_mat), nrow(.test_var_info(res)))
  expect_equal(ncol(expr_mat), nrow(.test_sample_info(res)))
})

# ----- Sample name extraction -----
test_that("sample names are correctly extracted from intensity columns", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  # Check that sample name matches the expected pattern
  expect_equal(.test_sample_info(res)$sample, "20241224-lxj-nglyco-h_1")
  expect_equal(colnames(.test_expr_mat(res)), "20241224-lxj-nglyco-h_1")
})

# ----- Error handling -----
test_that("TMT quantification throws error", {
  expect_error(
    read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "TMT"
    ),
    "TMT quantification is not supported yet"
  )
})

test_that("it handles O-linked glycan type", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      glycan_type = "O-GalNAc"
    )
  )
  expect_equal(.test_metadata(res)$glycan_type, "O-GalNAc")
})

# ----- Data consistency -----
test_that("all output data types are consistent", {
  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free"
    )
  )

  # Check data types of variable information
  expect_type(.test_var_info(res)$peptide, "character")
  expect_type(.test_var_info(res)$protein, "character")
  expect_type(.test_var_info(res)$protein_site, "integer")
  expect_type(.test_var_info(res)$peptide_site, "integer")

  # Check that all rows have data
  expect_true(all(nchar(.test_var_info(res)$peptide) > 0))
  expect_true(all(nchar(.test_var_info(res)$protein) > 0))
  expect_true(all(!is.na(.test_var_info(res)$protein_site)))
  expect_true(all(!is.na(.test_var_info(res)$peptide_site)))
  expect_true(all(.test_var_info(res)$protein_site > 0))
  expect_true(all(.test_var_info(res)$peptide_site > 0))
})

# ----- Sample name converter -----
test_that("sample name converter works", {
  sample_name_converter <- function(x) {
    stringr::str_replace(x, "20241224-lxj-nglyco-h_1", "Sample_1")
  }

  suppressMessages(
    res <- read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      sample_name_converter = sample_name_converter
    )
  )

  expect_equal(colnames(.test_expr_mat(res)), "Sample_1")
  expect_equal(.test_sample_info(res)$sample, "Sample_1")
})

# ----- Parameter validation (minimal, since common validation is tested in pglyco3) -----
test_that("it validates quant_method parameter", {
  expect_error(
    read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "invalid_method"
    ),
    " Must be element of set"
  )
})

test_that("it validates glycan_type parameter", {
  expect_error(
    read_byonic_pglycoquant(
      test_path("data/byonic-pglycoquant-LFQ-result.list"),
      quant_method = "label-free",
      glycan_type = "invalid_type"
    ),
    " Must be element of set"
  )
})

# ----- PSM aggregation tests -----
test_that("PSMs are correctly aggregated to glycopeptides", {
  # Read data and double the rows to simulate multiple PSMs per glycopeptide
  suppressMessages(
    res <- read_byonic_pglycoquant(test_path(
      "data/byonic-pglycoquant-LFQ-result.list"
    ))
  )
  # Check that aggregation worked correctly
  expect_equal(nrow(.test_var_info(res)), 5)
})
