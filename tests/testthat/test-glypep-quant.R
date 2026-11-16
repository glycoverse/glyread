test_that("site exports use the same row and sample annotations as Glyco-Decipher", {
  fp <- test_path("data/glyco-decipher-result.csv")
  site <- read_glypep_quant(fp)
  expect_equal(site, read_glyco_decipher(fp))
})

test_that("glycopeptide exports preserve modifications, candidates and abundances", {
  fp <- withr::local_tempfile(fileext = ".csv")
  writeLines(
    c(
      "Peptide,Glycan,Protein,HeLa_01,HeLa_02",
      '"ACVAGFTNITFHPQNK 2,Carbamidomethyl(+57.0215);",Hex(8)HexNAc(2),P1@135;P2@22;,10,',
      "ACVAGFTNITFHPQNK,Hex(8)HexNAc(2)+179.08,P1@135;,0,20",
      "ACVAGFTNITFHPQNK,Hex(8)HexNAc(2),P1@135;P1@142;,12,13",
      "EEQFNSTFR,Hex(3)HexNAc(3)Fuc(1),P01859,883.8,136.3"
    ),
    fp
  )
  x <- read_glypep_quant(
    fp,
    sample_info = data.frame(sample = c("02", "01")),
    sample_name_converter = \(s) sub("HeLa_", "", s)
  )
  expect_s4_class(x, "GlycoproteomicSE")
  expect_identical(colnames(x), c("02", "01"))
  expect_equal(
    unname(.test_expr_mat(x)),
    matrix(c(NA, 20, 13, 136.3, 10, NA, 12, 883.8), ncol = 2)
  )
  info <- .test_var_info(x)
  expect_identical(info$peptide, c(rep("ACVAGFTNITFHPQNK", 3), "EEQFNSTFR"))
  expect_identical(
    info$reported_peptide[1],
    "ACVAGFTNITFHPQNK 2,Carbamidomethyl(+57.0215);"
  )
  expect_identical(info$protein, c("P1", "P1", "P1", "P01859"))
  expect_identical(info$protein_group, c("P1;P2", "P1", "P1", "P01859"))
  expect_identical(info$protein_site, c(135L, 135L, NA_integer_, NA_integer_))
  expect_identical(info$reported_protein[1], "P1@135;P2@22;")
  expect_identical(info$reported_glycan[2], "Hex(8)HexNAc(2)+179.08")
  expect_equal(info$glycan_composition[1], info$glycan_composition[2])
  expect_identical(anyDuplicated(rownames(x)), 0L)
})

test_that("GlyPep-Quant requires a recognized header and at least one sample", {
  fp <- withr::local_tempfile(fileext = ".csv")
  writeLines("Peptide,Glycan,A", fp)
  expect_error(read_glypep_quant(fp), "Peptide, Glycan and Protein")
  writeLines("Protein,Peptide,Glycan", fp)
  expect_error(read_glypep_quant(fp), "sample abundance column")
})
