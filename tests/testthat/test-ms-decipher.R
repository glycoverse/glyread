.ms_decipher_fixture <- function() {
  test_path("data/ms-decipher-O_Search_15glycans.xlsx")
}

test_that("MS-Decipher imports the unchanged author workbook as binary detection", {
  raw <- .read_ms_decipher_df(.ms_decipher_fixture())
  expect_equal(nrow(raw), 2642L)
  expect_equal(
    readxl::excel_sheets(.ms_decipher_fixture()),
    c("Parameter", "PSM", "Peptide", "Protein")
  )
  result <- suppressMessages(read_ms_decipher(.ms_decipher_fixture()))
  expect_s4_class(result, "GlycoproteomicSE")
  expect_equal(dim(result), c(222L, 3L))
  mat <- .test_expr_mat(result)
  expect_true(all(mat %in% c(0, 1)))
  expect_equal(unname(colSums(mat)), c(190, 168, 174))
  expect_equal(colnames(mat), paste0("lly_QE_20201001_hplcogp_", 1:3))
  info <- .test_var_info(result)
  expect_true(all(is.na(info$peptide_site)))
  expect_true(all(is.na(info$protein_site)))
  expect_type(info$protein_site, "integer")
  expect_s3_class(info$glycan_composition, "glyrepr_composition")
  expect_equal(
    as.character(info$glycan_composition[1]),
    "Hex(2)HexNAc(2)NeuAc(2)"
  )
  expect_equal(info$peptide[1], "IEETTMTTQTPAPIQAPSAILPLPGQSVER")
  expect_equal(unname(mat[1, ]), c(1, 1, 1))
  expect_true("P00738;P00739" %in% info$protein_group)
  expect_false(any(c("gene", "glycan_structure") %in% names(info)))
  .expect_reader_metadata(
    result,
    list(
      exp_type = "glycoproteomics",
      glycan_type = "O-GalNAc",
      quant_method = "label-free"
    )
  )
})

test_that("MS-Decipher protein inference is independent of repeated PSMs", {
  raw <- .read_ms_decipher_df(.ms_decipher_fixture())
  original <- .tidy_ms_decipher(raw)
  repeated <- .tidy_ms_decipher(rbind(raw, raw[rep(1, 50), ]))
  expect_equal(dplyr::distinct(repeated), dplyr::distinct(original))
  groups <- original[original$protein_group == "P00738;P00739", ]
  expect_true(all(groups$protein %in% c("P00738", "P00739")))
})

test_that("MS-Decipher sample names and metadata can be supplied", {
  converter <- \(x) sub("lly_QE_20201001_hplcogp_", "sample", x)
  info <- data.frame(sample = paste0("sample", 3:1), group = "serum")
  result <- suppressMessages(read_ms_decipher(
    .ms_decipher_fixture(),
    sample_info = info,
    sample_name_converter = converter
  ))
  expect_setequal(colnames(result), paste0("sample", 1:3))
  expect_true(all(.test_sample_info(result)$group == "serum"))
})

test_that("MS-Decipher excludes decoys and rejects ambiguous or malformed annotations", {
  raw <- .read_ms_decipher_df(.ms_decipher_fixture())
  # Mutations of real source rows exercise validation, not alternate original files.
  decoy <- raw[1, ]
  decoy[["Hit Type"]] <- "DECOY"
  expect_equal(.tidy_ms_decipher(rbind(raw, decoy)), .tidy_ms_decipher(raw))
  expect_error(.tidy_ms_decipher(decoy), "no target")
  expect_error(.tidy_ms_decipher(raw[, -1]), "missing columns")
  bad <- raw[1, ]
  bad[["Delta Name"]] <- "unresolved delta"
  expect_error(.tidy_ms_decipher(bad), "composition")
  bad <- raw[1, ]
  bad$Proteins <- "P00738/P00739"
  expect_error(.tidy_ms_decipher(bad), "protein group")
  bad$Proteins <- "P00738::"
  expect_error(.tidy_ms_decipher(bad), "protein group")
  bad <- raw[1, ]
  bad[["Hit Type"]] <- "UNKNOWN"
  expect_error(.tidy_ms_decipher(bad), "TARGET or DECOY")
  bad <- raw[1, ]
  bad[["MS File"]] <- NA_character_
  expect_error(.tidy_ms_decipher(bad), "MS File")
})

test_that("MS-Decipher keeps decoy-only samples as zero detection columns", {
  raw <- .read_ms_decipher_df(.ms_decipher_fixture())
  decoy <- raw[1, ]
  decoy[["Hit Type"]] <- "DECOY"
  decoy[["MS File"]] <- "decoy_only_sample"
  testthat::local_mocked_bindings(.read_ms_decipher_df = function(fp) {
    rbind(raw, decoy)
  })
  result <- suppressMessages(read_ms_decipher(.ms_decipher_fixture()))
  expect_equal(dim(result), c(222L, 4L))
  expect_true(all(.test_expr_mat(result)[, "decoy_only_sample"] == 0))
})
