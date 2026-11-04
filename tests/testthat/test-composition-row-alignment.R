test_that("composition parsers retain missing rows when restoring duplicates", {
  parsers <- list(
    pglyco = list(parse = .convert_pglyco3_comp, valid = c("H(1)", "N(2)")),
    finder = list(parse = .parse_glycan_finder_composition,
                  valid = c("(Hex)1", "(HexNAc)2"))
  )
  expected <- glyrepr::as_glycan_composition(c(
    "Hex(1)", NA, "HexNAc(2)", NA, "Hex(1)", NA, "HexNAc(2)"
  ))

  for (parser in parsers) {
    x <- c(parser$valid[1], NA_character_, parser$valid[2], "",
           parser$valid[1], NA_character_, parser$valid[2])
    result <- parser$parse(x)
    expect_length(result, length(x))
    expect_equal(result, expected)
    expect_equal(parser$parse(c(NA_character_, NA_character_)),
                 glyrepr::glycan_composition(NA, NA))
    expect_equal(parser$parse(character()), glyrepr::glycan_composition())
  }
})
