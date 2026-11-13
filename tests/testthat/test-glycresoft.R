glycresoft_fixture <- function(name) {
  test_path("data", "glycresoft", paste0(name, ".csv"))
}
glycresoft_read <- function(names = c("HCC_2", "HCC_3", "Normal_3")) {
  read_glycresoft(
    vapply(names, glycresoft_fixture, ""),
    sample_name_converter = function(x) sub(" @.*$", "", x)
  )
}
glycresoft_mutate <- function(name, fn) {
  data <- readr::read_csv(
    glycresoft_fixture(name),
    col_types = readr::cols(.default = "c")
  )
  path <- tempfile(fileext = ".csv")
  readr::write_csv(fn(data), path)
  path
}

test_that("real shared groups retain one abundance and correct one-based sites", {
  x <- glycresoft_read("HCC_2")
  expect_s4_class(x, "GlycoproteomicSE")
  expect_true(methods::validObject(x))
  expect_equal(nrow(x), 4)
  info <- .test_var_info(x)
  abundance <- .test_expr_mat(x)[, 1]
  index <- which(info$peptide == "EEQFNSTFR")
  expect_equal(info$protein_site[index], 176L)
  expect_equal(info$peptide_site[index], 5L)
  expect_equal(unname(abundance[index]), 1426909927.8125)
  expect_equal(
    sum(abundance, na.rm = TRUE),
    sum(c(
      1426909927.8125,
      6478411356.578125,
      2036066037.125,
      as.numeric(readr::read_csv(
        glycresoft_fixture("HCC_2"),
        show_col_types = FALSE
      )$total_signal[
        readr::read_csv(
          glycresoft_fixture("HCC_2"),
          show_col_types = FALSE
        )$group_id ==
          118
      ][1])
    ))
  )
  ambiguous <- which(vapply(
    info$source_annotations,
    function(z) any(z$group_id == "81"),
    TRUE
  ))
  expect_true(is.na(info$protein_site[ambiguous]))
  expect_equal(nrow(info$source_annotations[[ambiguous]]), 4)
  expect_match(info$protein_group[ambiguous], "P00738")
  expect_match(info$protein_group[ambiguous], "P00739")
  expect_equal(.test_metadata(x)$abundance_source, "GlycReSoft total_signal")
})

test_that("different unknown-composition candidate sets remain separate and align", {
  x <- glycresoft_read()
  info <- .test_var_info(x)
  index <- which(is.na(info$glycan_composition))
  expect_length(index, 2)
  expect_true(all(lengths(info$candidate_compositions[index]) == 2))
  mat <- .test_expr_mat(x)
  normal <- grep("Normal_3", colnames(mat))
  hcc2 <- grep("HCC_2", colnames(mat))
  hcc3 <- grep("HCC_3", colnames(mat))
  one <- index[which(mat[index, normal] == 554962798.75)]
  two <- index[which(mat[index, normal] == 1998217504.125)]
  expect_length(one, 1)
  expect_length(two, 1)
  expect_false(identical(one, two))
  expect_false(is.na(mat[one, hcc3]))
  expect_true(is.na(mat[one, hcc2]))
  expect_equal(mat[two, hcc2], 2036066037.125)
  expect_false(identical(
    info$feature_signature[one],
    info$feature_signature[two]
  ))
})

test_that("real yeast output keeps Neu5Gc and zero-abundance identities", {
  x <- glycresoft_read("cwq_mix2-1")
  info <- .test_var_info(x)
  first <- which(info$peptide == "LGNNLTR")
  expect_equal(info$protein_site[first], 209L)
  expect_equal(info$protein[first], "Q9C0Y4")
  expect_equal(info$peptide_site[first], 4L)
  expect_true(any(grepl("NeuGc", as.character(info$glycan_composition))))
  expect_true(any(is.na(.test_expr_mat(x))))
})

test_that("sample-local groups are not global feature identifiers", {
  path <- glycresoft_mutate("HCC_2", function(d) {
    d <- d[d$group_id == "1", ]
    second <- d
    second$group_id <- "9001"
    rbind(d, second)
  })
  x <- read_glycresoft(path)
  expect_equal(nrow(x), 1)
  expect_equal(as.numeric(.test_expr_mat(x)), 2 * 1426909927.8125)
  expect_equal(nrow(.test_var_info(x)$source_annotations[[1]]), 2)
})

test_that("sample naming, ordering and duplicate analysis validation work", {
  x <- glycresoft_read()
  samples <- rev(colnames(x))
  y <- read_glycresoft(
    vapply(c("HCC_2", "HCC_3", "Normal_3"), glycresoft_fixture, ""),
    sample_info = data.frame(sample = samples),
    sample_name_converter = function(x) sub(" @.*$", "", x)
  )
  expect_identical(colnames(y), samples)
  expect_equal(.test_expr_mat(y), .test_expr_mat(x)[, samples])
  path <- glycresoft_mutate("HCC_2", identity)
  expect_error(
    read_glycresoft(c(glycresoft_fixture("HCC_2"), path)),
    "multiple files"
  )
  expect_error(read_glycresoft(
    vapply(c("HCC_2", "HCC_3"), glycresoft_fixture, ""),
    sample_name_converter = function(x) rep("same", length(x))
  ))
})

test_that("invalid or unsupported exports fail without guessing shared groups", {
  expect_error(
    read_glycresoft(glycresoft_mutate("HCC_2", function(x) {
      x$group_id <- "-1"
      x
    })),
    "group_id"
  )
  expect_error(
    read_glycresoft(glycresoft_mutate("HCC_2", function(x) {
      x$n_glycosylation_sites <- "175:NST"
      x
    })),
    "Unsupported n_glycosylation_sites"
  )
  expect_error(
    read_glycresoft(glycresoft_mutate("HCC_2", function(x) {
      x$total_signal[x$group_id == "81"][1] <- "-1"
      x
    })),
    "total_signal"
  )
  expect_error(
    read_glycresoft(glycresoft_mutate("HCC_2", function(x) {
      x$total_signal[which(x$group_id == "81")[1]] <- "100"
      x
    })),
    "Conflicting total_signal"
  )
  expect_error(
    read_glycresoft(glycresoft_mutate("HCC_2", function(x) {
      x$glycopeptide <- sub("Fuc:", "Unknown:", x$glycopeptide)
      x
    })),
    "Unsupported GlycReSoft monosaccharide"
  )
  expect_error(
    read_glycresoft(glycresoft_mutate("HCC_2", function(x) {
      x$group_id <- NULL
      x
    })),
    "identification CSV"
  )
})

test_that("less-than-certain single candidates retain unknown sites", {
  path <- glycresoft_mutate("HCC_2", function(x) {
    x <- x[x$group_id == "1", ]
    x$localizations <- "N-Glycosylation:4:0.99"
    x
  })
  x <- read_glycresoft(path)
  expect_true(is.na(.test_var_info(x)$protein_site))
  expect_true(is.na(.test_var_info(x)$peptide_site))
  expect_equal(as.numeric(.test_expr_mat(x)), 1426909927.8125)
})

test_that("parsimony aligns subset protein mappings without losing paired sites", {
  original <- readr::read_csv(
    glycresoft_fixture("cwq_mix2-1"),
    col_types = readr::cols(.default = "c")
  )
  original <- original[original$group_id == "1", ]
  first <- tempfile(fileext = ".csv")
  second <- tempfile(fileext = ".csv")
  readr::write_csv(original, first)
  subset <- original[grepl("Q9C0Y4", original$protein_name), ]
  subset$analysis <- "another.analysis"
  subset$group_id <- "500"
  readr::write_csv(subset, second)
  x <- read_glycresoft(c(first, second))
  expect_equal(nrow(x), 1)
  expect_equal(.test_var_info(x)$protein, "Q9C0Y4")
  expect_equal(.test_var_info(x)$protein_site, 209L)
  expect_equal(unname(.test_expr_mat(x)[1, ]), rep(106946901.2734375, 2))
})

test_that("different analysis and group tuples cannot collide", {
  path <- glycresoft_mutate("HCC_2", function(x) {
    x <- x[x$group_id == "1", ]
    y <- x
    x$analysis <- "a.b"
    y$analysis <- "a"
    y$group_id <- "2"
    rbind(x, y)
  })
  x <- read_glycresoft(path)
  expect_equal(dim(x), c(1L, 2L))
  expect_equal(unname(.test_expr_mat(x)[1, ]), rep(1426909927.8125, 2))
})

test_that("multioccupied sites stay unknown rather than splitting a glycan", {
  path <- glycresoft_mutate("HCC_2", function(x) {
    x <- x[x$group_id == "118", ][1, ]
    x$glycopeptide <- "N(N-Glycosylation)VTN(N-Glycosylation)AT{Hex:5; HexNAc:4}"
    x$peptide_start <- "10"
    x$peptide_end <- "16"
    x$n_glycosylation_sites <- "10;13"
    x$localizations <- "N-Glycosylation:0:1;N-Glycosylation:3:1"
    x
  })
  x <- read_glycresoft(path)
  expect_true(is.na(.test_var_info(x)$protein_site))
  expect_true(is.na(.test_var_info(x)$peptide_site))
  expect_equal(nrow(x), 1)
})

test_that("glycan naming matches the other glyread readers", {
  parsed <- .glycresoft_parse(c(
    "N(N-Glycosylation)ST{Fuc:1; Hex:3; HexNAc:4; Neu5Ac:2}",
    "N(N-Glycosylation)ST{Hex:3; HexNAc:4; Neu5Gc:1}"
  ))
  expect_identical(
    parsed$compositions,
    as.character(glyrepr::as_glycan_composition(c(
      "Hex(3)HexNAc(4)dHex(1)NeuAc(2)",
      "Hex(3)HexNAc(4)NeuGc(1)"
    )))
  )
})

test_that("uncertain localization in one sample is not hidden by another sample", {
  first <- glycresoft_mutate("HCC_2", function(x) x[x$group_id == "1", ])
  second <- glycresoft_mutate("HCC_2", function(x) {
    x <- x[x$group_id == "1", ]
    x$analysis <- "less.certain"
    x$localizations <- "N-Glycosylation:4:0.99"
    x
  })
  x <- read_glycresoft(c(first, second))
  expect_equal(nrow(x), 1)
  expect_true(is.na(.test_var_info(x)$protein_site))
  expect_equal(unname(.test_expr_mat(x)[1, ]), rep(1426909927.8125, 2))
})

test_that("semicolons in FASTA descriptions do not change protein accessions", {
  path <- glycresoft_mutate("HCC_2", function(x) {
    x$protein_name <- paste0(x$protein_name, "; additional description")
    x
  })
  original <- read_glycresoft(glycresoft_fixture("HCC_2"))
  changed <- read_glycresoft(path)
  expect_identical(
    .test_var_info(original)$protein,
    .test_var_info(changed)$protein
  )
  expect_identical(
    .test_var_info(original)$protein_site,
    .test_var_info(changed)$protein_site
  )
})
