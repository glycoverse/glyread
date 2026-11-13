.gprodia_fixture <- function(name = "serum_uis") {
  test_path(paste0("data/gprodia-", name, "_aligned.tsv"))
}

.gprodia_raw <- function(name = "serum_uis") {
  readr::read_tsv(
    .gprodia_fixture(name),
    col_types = readr::cols(.default = "c")
  )
}

test_that("GproDIA reads original author results and preserves MS2 totals", {
  for (name in c("serum_uis", "fissionyeast")) {
    raw <- .gprodia_raw(name)
    result <- suppressMessages(read_gprodia(.gprodia_fixture(name)))
    expect_s4_class(result, "GlycoproteomicSE")
    info <- .test_var_info(result)
    expect_s3_class(info$glycan_composition, "glyrepr_composition")
    expect_false(any(grepl("J", info$peptide)))
    expect_false("glycan_structure" %in% names(info))
    expect_setequal(colnames(result), unique(raw$filename))
    target <- raw$decoy == "0" &
      raw$decoy_peptide == "0" &
      raw$decoy_glycan == "0"
    expected <- tapply(
      as.numeric(raw$Intensity[target]),
      raw$filename[target],
      sum
    )
    actual <- colSums(.test_expr_mat(result), na.rm = TRUE)
    expect_equal(unname(actual[names(expected)]), as.numeric(expected))
    .expect_reader_metadata(
      result,
      list(
        exp_type = "glycoproteomics",
        glycan_type = "N",
        quant_method = "label-free"
      )
    )
  }
})

test_that("GproDIA resolves source slash-separated proteins and paired sites", {
  tidy <- suppressMessages(.tidy_gprodia(.gprodia_raw()))
  rows <- tidy[tidy$protein_group == "P01871;P0DOX6", ]
  expect_gt(nrow(rows), 0)
  expect_true(all(rows$protein == "P01871"))
  expect_true(all(rows$protein_site == 209L))
  # Semicolon is also emitted by the upstream SQL GROUP_CONCAT exporter.
  raw <- .gprodia_raw()
  raw$ProteinName <- gsub("/", ";", raw$ProteinName, fixed = TRUE)
  raw$ProteinGlycoSite <- gsub("/", ";", raw$ProteinGlycoSite, fixed = TRUE)
  expect_equal(suppressMessages(.tidy_gprodia(raw)), tidy)
})

test_that("GproDIA selects the lowest-score peak group without double counting", {
  raw <- .gprodia_raw()
  extra <- raw[1, ]
  extra$m_score <- "0.99"
  extra$Intensity <- "999999999"
  expect_equal(
    suppressMessages(.tidy_gprodia(rbind(raw, extra))),
    suppressMessages(.tidy_gprodia(raw))
  )
  # Real serum output includes slightly negative floating-point q-values.
  expect_true(any(as.numeric(raw$m_score) < 0))
})

test_that("GproDIA accepts sample information and sample name conversion", {
  convert <- \(x) sub("[.]mzML$", "", x)
  samples <- convert(unique(.gprodia_raw()$filename))
  info <- data.frame(sample = rev(samples), group = "serum")
  result <- suppressMessages(read_gprodia(
    .gprodia_fixture(),
    sample_info = info,
    sample_name_converter = convert
  ))
  expect_setequal(colnames(result), samples)
  expect_true(all(.test_sample_info(result)$group == "serum"))
  expect_error(
    read_gprodia(.gprodia_fixture(), quant_method = "TMT"),
    "not supported"
  )
})

test_that("GproDIA rejects malformed reports rather than guessing", {
  raw <- .gprodia_raw()
  expect_error(.tidy_gprodia(raw[, -1]), "missing columns")
  bad <- raw
  bad$ProteinGlycoSite[1] <- "1/2"
  expect_error(.tidy_gprodia(bad), "paired by position")
  bad <- raw
  bad$GlycanComposition[1] <- "X(2)"
  expect_error(.tidy_gprodia(bad), "composition notation")
  bad <- raw
  bad$Intensity[1] <- "invalid"
  expect_error(.tidy_gprodia(bad), "invalid numeric")
  bad <- raw
  bad$GlycanSite[1] <- "0"
  expect_error(.tidy_gprodia(bad), "one-based")
  bad <- raw
  bad$decoy[1] <- NA_character_
  expect_error(.tidy_gprodia(bad), "only 0 or 1")
  bad <- raw
  bad$decoy <- "1"
  expect_error(.tidy_gprodia(bad), "no target rows")
})


test_that("GproDIA sums real modified and unmodified precursors by composition", {
  result <- suppressMessages(read_gprodia(.gprodia_fixture()))
  info <- .test_var_info(result)
  index <- which(info$peptide == "GLTFQQNASSMCVPDQDTAIR")
  expect_length(index, 1)
  expect_equal(info$peptide_site[index], 7L)
  expect_equal(info$protein_site[index], 209L)
  expect_equal(
    unname(.test_expr_mat(result)[
      index,
      "20200615_serum_NGlycP_1hDIA_rep1.mzML"
    ]),
    22237400 + 1439100
  )
})


test_that("GproDIA refuses run identifiers reused across different files", {
  raw <- .gprodia_raw()
  raw$run_id <- raw$run_id[1]
  expect_error(.tidy_gprodia(raw), "exactly one filename")
})
