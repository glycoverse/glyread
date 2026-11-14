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
  .expect_reader_metadata(
    res,
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
      Composition = "HexNAc(4)Hex(5)Fuc(1)NeuAc(1),HexNAc(1)Fuc(1)",
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
  expect_equal(var_info$protein_site, c(243L, 247L))
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

test_that("mass matches locate glycosylation on C, S and T among other modifications", {
  path <- withr::local_tempfile(fileext = ".list")
  readr::write_tsv(
    tibble::tibble(
      Peptide = c(
        "K.N[+0.98402]C[+203.07937]AR.G",
        "R.M[+15.99492]AST[+162.05282]S[+203.07937]R.K"
      ),
      Composition = c("HexNAc(1)", "HexNAc(1),Hex(1)"),
      Position = c(10L, 20L),
      `Protein Name` = ">sp|P12345|TEST",
      `Intensity(Sample1)` = c(10, 20)
    ),
    path
  )
  res <- read_byonic_pglycoquant(path, orgdb = "missing.OrgDb")
  info <- dplyr::arrange(.test_var_info(res), .data$protein_site)
  expect_equal(info$peptide_site, c(2L, 4L, 5L))
  expect_equal(info$protein_site, c(11L, 23L, 24L))
  expect_equal(
    as.character(info$glycan_composition),
    c("HexNAc(1)", "Hex(1)", "HexNAc(1)")
  )
})

test_that("partial and ambiguous matches keep known fields and warn once per import", {
  path <- withr::local_tempfile(fileext = ".list")
  readr::write_tsv(
    tibble::tibble(
      Peptide = c(
        "K.S[+203.07937]T[+203.07937]C[+162.05282].R",
        "R.S[+162.05282]T[+162.05282]N[+203.07937].K",
        "R.N[+0.98402]AST.K",
        "R.S[+203.07937]T[+203.07937].K"
      ),
      Composition = c(
        "Hex(1),HexNAc(1)",
        "Gal(1),Glc(1),HexNAc(1)",
        "HexNAc(1),HexNAc(1)",
        "HexNAc(1),HexNAc(1)"
      ),
      Position = c(10L, 20L, 30L, 40L),
      `Protein Name` = ">sp|P12345|TEST",
      `Intensity(Sample1)` = c(100, 200, 300, 400)
    ),
    path
  )
  expect_warning(
    res <- read_byonic_pglycoquant(path, orgdb = "missing.OrgDb"),
    "3 of 4",
    class = "glyread_unresolved_byonic_sites"
  )
  info <- .test_var_info(res)
  values <- as.numeric(.test_expr_mat(res)[, "Sample1"])
  info <- info[order(values, info$peptide_site), ]
  expect_equal(sort(values), c(100, 100, 200, 200, 200, 300, 400, 400))
  expect_equal(info$peptide_site, c(3L, NA, 1L, 2L, 3L, NA, 1L, 2L))
  expect_equal(info$protein_site, c(12L, NA, 20L, 21L, 22L, NA, 40L, 41L))
  expect_equal(
    is.na(info$glycan_composition),
    c(FALSE, FALSE, TRUE, TRUE, rep(FALSE, 4))
  )
  expect_equal(
    as.character(info$glycan_composition[!is.na(info$glycan_composition)]),
    c("Hex(1)", rep("HexNAc(1)", 5))
  )
})

test_that("competing glycans use each site once and preserve an independent unmatched glycan", {
  df <- tibble::tibble(
    Peptide = "K.S[+57.02146][+162.05282]T[+203.08050].R",
    Composition = "Gal(1),Glc(1),HexNAc(1)",
    Position = 8L
  )
  expect_warning(
    pairs <- .expand_byonic_pglycoquant_multisite_rows(df),
    class = "glyread_unresolved_byonic_sites"
  )
  # The sole 162-Da site is occupied. Either hexose could occupy it.
  expect_equal(pairs$protein_site, c(8L, NA, NA))
  expect_true(is.na(pairs$glycan_composition[1]))
  # The other hexose has unknown identity and site. HexNAc is outside tolerance.
  expect_setequal(pairs$glycan_composition[-1], c(NA_character_, "HexNAc(1)"))
})

test_that("unmarked and unsupported assignments do not create sites", {
  df <- tibble::tibble(
    Peptide = c("K.nGTR.G", "K.N[+203.07937]T.R", "K.AT.R"),
    Composition = c("HexNAc(1)", "Unknown(1)", NA_character_),
    Position = 10L
  )
  expect_warning(
    pairs <- .handle_byonic_pglycoquant_multisite_rows(df, "drop"),
    "3 of 3",
    class = "glyread_unresolved_byonic_sites"
  )
  expect_equal(pairs$peptide_site, rep(NA_integer_, 3))
  expect_equal(pairs$glycan_composition, c("HexNAc(1)", NA, NA))
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

test_that("site tables import reported sites and sample abundances directly", {
  path <- withr::local_tempfile(fileext = ".list")
  writeLines(
    c(
      "Protein\tSite\tGlycan(N,H,A,G,F)\tEmpty_Separator\trun-1\trun-2\t",
      ">sp|P19652-1|A1AG2_HUMAN\t  102\t[2  3  1  0  1]\t\t0.0\t32.5\t",
      ">sp|P19652-1|A1AG2_HUMAN\t  102\t[2  3  0  1  0]\t\t11.25\t\t"
    ),
    path
  )

  res <- read_byonic_pglycoquant(
    path,
    sample_info = data.frame(sample = c("2", "1")),
    sample_name_converter = function(x) sub("run-", "", x),
    orgdb = "missing.OrgDb",
    multisite = "drop"
  )
  info <- .test_var_info(res)
  expect_equal(info$protein, rep("P19652-1", 2))
  expect_identical(info$protein_site, c(102L, 102L))
  expect_false(any(c("peptide", "peptide_site") %in% names(info)))
  expect_s3_class(info$glycan_composition, "glyrepr_composition")
  expect_equal(
    as.character(info$glycan_composition),
    c(
      "Hex(3)HexNAc(2)dHex(1)NeuAc(1)",
      "Hex(3)HexNAc(2)NeuGc(1)"
    )
  )
  mat <- .test_expr_mat(res)
  expect_identical(colnames(mat), c("2", "1"))
  expect_equal(unname(mat), matrix(c(32.5, NA, NA, 11.25), nrow = 2))
})

test_that("site compositions require a count for every declared sugar", {
  expect_error(
    .parse_byonic_site_composition("[3  2]", "Glycan(H,N,F)"),
    "one nonnegative count per monosaccharide"
  )
})
