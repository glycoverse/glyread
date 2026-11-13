.pp_fixture <- function() {
  system.file("extdata", "protein-prospector.tsv", package = "glyread")
}

test_that("real Protein Prospector records retain precursor intensity semantics", {
  df <- .read_protein_prospector_df(.pp_fixture())
  tidy <- .tidy_protein_prospector(df)
  res <- suppressMessages(read_protein_prospector(
    .pp_fixture(),
    orgdb = "not.installed"
  ))
  expect_s4_class(res, "GlycoproteomicSE")
  expect_equal(.test_metadata(res)$abundance_source, "PSM precursor intensity")
  expect_equal(.test_metadata(res)$quant_method, "label-free")
  expect_equal(
    sum(.test_expr_mat(res), na.rm = TRUE),
    sum(as.numeric(df$Intensity))
  )
  expect_equal(tidy$protein_site[1], 509L)
  expect_equal(tidy$peptide_site[1], 1L)
  expect_equal(tidy$value[1], 88202)
  expect_true(anyNA(tidy$peptide_site))
  expect_true(any(grepl("dHex", tidy$glycan_composition)))
  expect_setequal(colnames(res), df$Fraction)
  expect_equal(
    nrow(res),
    nrow(dplyr::distinct(
      tidy,
      peptide,
      peptide_site,
      protein,
      protein_site,
      glycan_composition
    ))
  )
  ambiguous <- grepl("[|]", df$Mods) | endsWith(df$Mods, "@")
  expect_true(all(is.na(tidy$peptide_site[ambiguous])))
  expect_true(all(is.na(tidy$protein_site[ambiguous])))
})

test_that("PSMs are summed per feature and fraction", {
  tidy <- .tidy_protein_prospector(.read_protein_prospector_df(.pp_fixture()))
  expected <- .aggregate_long(tidy)
  res <- suppressMessages(read_protein_prospector(
    .pp_fixture(),
    orgdb = "not.installed"
  ))
  vars <- .test_var_info(res)
  mat <- .test_expr_mat(res)
  for (i in seq_len(nrow(expected))) {
    same_site <- (is.na(vars$peptide_site) & is.na(expected$peptide_site[i])) |
      (!is.na(vars$peptide_site) &
        !is.na(expected$peptide_site[i]) &
        vars$peptide_site == expected$peptide_site[i])
    match <- which(
      vars$peptide == expected$peptide[i] &
        vars$protein == expected$protein[i] &
        same_site &
        as.character(vars$glycan_composition) ==
          as.character(glyrepr::as_glycan_composition(expected$glycan_composition[
            i
          ]))
    )
    expect_length(match, 1)
    expect_equal(unname(mat[match, expected$sample[i]]), expected$value[i])
  }
})

test_that("sample metadata and sample names follow the reader interface", {
  res <- suppressMessages(read_protein_prospector(
    .pp_fixture(),
    orgdb = "not.installed",
    sample_name_converter = toupper
  ))
  expect_true(all(colnames(res) == toupper(colnames(res))))
  info <- data.frame(sample = rev(colnames(res)), batch = "A")
  other <- suppressMessages(read_protein_prospector(
    .pp_fixture(),
    sample_info = info,
    sample_name_converter = toupper,
    orgdb = "not.installed"
  ))
  expect_equal(colnames(other), info$sample)
  expect_equal(
    as.character(SummarizedExperiment::colData(other)$batch),
    info$batch
  )
})

test_that("invalid or unsupported fields are rejected rather than guessed", {
  df <- .read_protein_prospector_df(.pp_fixture())[1, ]
  bad <- df
  bad$Intensity <- "broken"
  expect_error(.tidy_protein_prospector(bad), "Intensity")
  bad <- df
  bad$Intensity <- "-1"
  expect_error(.tidy_protein_prospector(bad), "Intensity")
  bad <- df
  bad[["Start AA"]] <- "1.5"
  expect_error(.tidy_protein_prospector(bad), "Start AA")
  bad <- df
  bad$Mods <- "HexNAc2Hex4*@999"
  expect_error(.tidy_protein_prospector(bad), "within")
  bad <- df
  bad$Mods <- "HexNAc2Hex4*@1?"
  expect_error(.tidy_protein_prospector(bad), "syntax")
  bad <- df
  bad[["Acc Num"]] <- "P12821;P12345"
  expect_error(.tidy_protein_prospector(bad), "protein accessions")
  expect_error(
    .protein_prospector_composition("HexNAc2Unknown1"),
    "Unsupported"
  )
  expect_error(.protein_prospector_composition("Hex2Hex3"), "Repeated")
  fp <- tempfile(fileext = ".tsv")
  readr::write_tsv(df[, -1], fp)
  expect_error(.read_protein_prospector_df(fp), "requires")
})

test_that("aligned protein groups use the chosen protein's coordinate", {
  df <- .read_protein_prospector_df(.pp_fixture())
  df <- df[grepl(",", df[["Acc Num"]]), ]
  tidy <- .tidy_protein_prospector(df)
  expect_equal(tidy$protein, "O75477")
  expect_equal(tidy$protein_site, 108L)
  expect_equal(tidy$protein_group, "O75477;O94905")
  bad <- df
  bad[["Start AA"]] <- "108"
  expect_error(.tidy_protein_prospector(bad), "aligned")
  bad <- df
  bad[["Start AA"]] <- NA_character_
  expect_true(all(is.na(.tidy_protein_prospector(bad)$protein_site)))
})

test_that("missing and zero intensities remain missing after aggregation", {
  df <- .read_protein_prospector_df(.pp_fixture())[1, ]
  df$Intensity <- NA_character_
  expect_true(is.na(.aggregate_long(.tidy_protein_prospector(df))$value))
  df$Intensity <- "0"
  expect_true(is.na(.aggregate_long(.tidy_protein_prospector(df))$value))
  expect_equal(
    .protein_prospector_composition("HexNAc2HexFuc"),
    "HexNAc(2)Hex(1)dHex(1)"
  )
})
