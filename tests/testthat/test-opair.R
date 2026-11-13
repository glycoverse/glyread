.opair_fixture <- function() {
  x <- readr::read_lines(test_path("data/opair/oglyco.psmtsv"))
  x[1] <- sub("\\t+$", "", x[1])
  readr::read_tsv(
    I(paste(x, collapse = "\n")),
    col_types = readr::cols(.default = "c")
  )
}

.opair_write <- function(x) {
  fp <- tempfile(fileext = ".psmtsv")
  readr::write_tsv(x, fp, na = "")
  withr::defer(unlink(fp), envir = parent.frame())
  fp
}

test_that("native O-Pair output gives complete-glycopeptide detections", {
  x <- suppressMessages(read_opair(test_path("data/opair/oglyco.psmtsv")))
  expect_s4_class(x, "GlycoproteomicSE")
  expect_equal(dim(x), c(7L, 4L))
  expect_true(all(.test_expr_mat(x) %in% c(0, 1)))
  expect_true(any(.test_expr_mat(x) == 0))
  v <- .test_var_info(x)
  expect_s3_class(v$glycan_composition, "glyrepr_composition")
  expect_equal(v$protein_site[v$protein == "P01326"], c(79L, 79L))
  expect_true(all(is.na(v$protein_site[v$protein != "P01326"])))
  expect_true("A2AS86;P55095" %in% v$protein_group)
  expect_equal(
    as.character(v$glycan_composition[v$protein == "P01326"]),
    rep("Hex(1)HexNAc(1)NeuAc(2)", 2)
  )
  expect_true(all(
    c(
      "full_sequence",
      "peptide_localization",
      "protein_localization",
      "potential_localizations",
      "site_probabilities"
    ) %in%
      names(v)
  ))
  .expect_reader_metadata(
    x,
    list(
      exp_type = "glycoproteomics",
      glycan_type = "O-GalNAc",
      quant_method = "label-free"
    )
  )
})

test_that("repeated PSMs and probability changes do not create counts or features", {
  df <- .opair_fixture()
  a <- suppressMessages(read_opair(.opair_write(df)))
  extra <- df[1, ]
  extra[[
    "Localized Glycans with Protein Site Specific Probability"
  ]] <- "[79,H1N1A2,0.900]"
  b <- suppressMessages(read_opair(.opair_write(dplyr::bind_rows(
    df,
    df,
    extra
  ))))
  expect_equal(dim(a), dim(b))
  expect_equal(unname(.test_expr_mat(a)), unname(.test_expr_mat(b)))
  v <- .test_var_info(b)
  changed <- grepl("0.900", v$protein_localization, fixed = TRUE)
  expect_true(any(changed))
  expect_true(all(is.na(v$protein_site[changed])))
  expect_true(all(grepl(" || ", v$protein_localization[changed], fixed = TRUE)))
})

test_that("multiple glycans retain total composition without assigning it to sites", {
  df <- .opair_fixture()[1, ]
  # Derived edge case: two site-specific glycans share a whole-peptide total.
  df[["Plausible Number Of Glycans"]] <- "2"
  df[[
    "Localized Glycans with Protein Site Specific Probability"
  ]] <- "[79,H1N1A1,1.000][80,A1,1.000]"
  x <- suppressMessages(read_opair(.opair_write(df)))
  expect_equal(nrow(x), 1L)
  expect_true(is.na(.test_var_info(x)$protein_site))
  expect_equal(
    as.character(.test_var_info(x)$glycan_composition),
    "Hex(1)HexNAc(1)NeuAc(2)"
  )
  expect_equal(
    .test_var_info(x)$protein_localization,
    df[["Localized Glycans with Protein Site Specific Probability"]]
  )
})

test_that("target and q filtering retain all input samples and respect sample order", {
  df <- .opair_fixture()
  extra <- df[1:2, ]
  extra[["File Name"]] <- c("decoy-only", "high-q-only")
  extra[["Decoy/Contaminant/Target"]][1] <- "D"
  extra$QValue[2] <- "0.5"
  df <- dplyr::bind_rows(df, extra)
  info <- data.frame(sample = rev(unique(df[["File Name"]])))
  x <- suppressMessages(read_opair(
    .opair_write(df),
    sample_info = info,
    max_q_value = 0.01
  ))
  expect_identical(colnames(x), info$sample)
  expect_true(all(.test_expr_mat(x)[, c("decoy-only", "high-q-only")] == 0))
  y <- suppressMessages(read_opair(.opair_write(df)))
  expect_true(any(.test_expr_mat(y)[, "high-q-only"] == 1))
  expect_error(
    read_opair(.opair_write(df), sample_info = data.frame(sample = "missing")),
    "must match"
  )
})

test_that("invalid native outputs fail clearly", {
  df <- .opair_fixture()
  expect_error(read_opair(.opair_write(df), max_q_value = -1))
  expect_error(read_opair(.opair_write(df[, names(df) != "QValue"])), "QValue")
  bad <- df
  bad$QValue[1] <- "not-a-number"
  expect_error(suppressMessages(read_opair(.opair_write(bad))), "QValue must")
  bad <- df
  bad[["Decoy/Contaminant/Target"]] <- "D"
  expect_error(suppressMessages(read_opair(.opair_write(bad))), "No target")
  bad <- df
  bad[["Plausible GlycanComposition"]][1] <- "H1N1CUSTOM1"
  expect_error(suppressMessages(read_opair(.opair_write(bad))), "Unsupported")
  bad <- df
  bad[["File Name"]][1] <- NA_character_
  expect_error(read_opair(.opair_write(bad)), "missing values")
  expect_error(read_opair(.opair_write(df[FALSE, ])), "no PSMs")
})

test_that("composition parser follows MetaMorpheus notation, not StrucGP S", {
  expect_equal(
    as.character(.parse_opair_comp("H5N4A2G1F1")),
    as.character(glyrepr::as_glycan_composition(
      "Hex(5)HexNAc(4)NeuAc(2)NeuGc(1)dHex(1)"
    ))
  )
  expect_error(.parse_opair_comp("H1N1S1"), "Unsupported")
  expect_error(.parse_opair_comp("H1H2"), "Repeated")
  expect_error(.parse_opair_comp(NA_character_), "missing")
})

test_that("O-Pair HTTP input matches the unchanged local fixture", {
  skip_if_not_installed("httpuv")
  root <- normalizePath(test_path("data/opair"))
  port <- tryCatch(httpuv::randomPort(), httpuv_unavailable_port = function(e) {
    skip("No local port is available")
  })
  server <- httpuv::startServer(
    "127.0.0.1",
    port,
    list(staticPaths = list("/" = root))
  )
  on.exit(server$stop(), add = TRUE)
  url <- paste0("http://127.0.0.1:", port, "/oglyco.psmtsv")
  expect_equal(
    suppressMessages(read_opair(url)),
    suppressMessages(read_opair(file.path(root, "oglyco.psmtsv")))
  )
})


test_that("leader inference precedes detection grouping across candidate groups", {
  df <- .opair_fixture()[rep(1, 2), ]
  df[["File Name"]] <- c("sample-one", "sample-two")
  df[["Protein Accession"]] <- c("P01326", "P01326|P55095")
  x <- suppressMessages(read_opair(.opair_write(df)))
  expect_equal(dim(x), c(1L, 2L))
  expect_equal(unname(.test_expr_mat(x)), matrix(1, 1, 2))
  expect_identical(.test_var_info(x)$protein, "P01326")
  expect_identical(.test_var_info(x)$protein_group, "P01326;P55095")
  expect_true(is.na(.test_var_info(x)$protein_site))
})


test_that("invalid classifications fail and empty decoy identity does not block targets", {
  df <- .opair_fixture()
  for (flag in c(NA_character_, "unexpected")) {
    bad <- df
    bad[["Decoy/Contaminant/Target"]][1] <- flag
    expect_error(read_opair(.opair_write(bad)), "only T, D or C")
  }
  decoy <- df[1, ]
  decoy[["Decoy/Contaminant/Target"]] <- "D"
  decoy[["Protein Accession"]] <- NA_character_
  decoy[["Base Sequence"]] <- NA_character_
  decoy[["Full Sequence"]] <- NA_character_
  expect_equal(
    suppressMessages(read_opair(.opair_write(dplyr::bind_rows(df, decoy)))),
    suppressMessages(read_opair(.opair_write(df)))
  )
})


test_that("repeating one spectrum cannot change leader inference", {
  df <- .opair_fixture()[1:3, ]
  df[["Protein Accession"]] <- c("P01326|P55095", "P01326", "P55095")
  a <- suppressMessages(read_opair(.opair_write(df)))
  b <- suppressMessages(read_opair(.opair_write(dplyr::bind_rows(
    df,
    df[rep(3, 20), ]
  ))))
  expect_equal(a, b)
})
