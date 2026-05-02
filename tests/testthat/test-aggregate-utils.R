test_that("aggregation ignores trace columns and sums glycopeptide values", {
  tidy_df <- tibble::tibble(
    peptide = "NTEST",
    peptide_site = 1L,
    protein = "P12345",
    protein_site = 10L,
    glycan_composition = "HexNAc(1)",
    gp_id = c("GP1", "GP2"),
    sample = "S1",
    value = c(10, 20)
  )

  res <- glyread:::.aggregate_long(tidy_df)

  expect_equal(nrow(res), 1L)
  expect_equal(res$value, 30)
  expect_false("gp_id" %in% colnames(res))
})

test_that("aggregation keeps sample-level values separate", {
  tidy_df <- tibble::tibble(
    peptide = "NTEST",
    peptide_site = 1L,
    protein = "P12345",
    protein_site = 10L,
    glycan_composition = "HexNAc(1)",
    gp_id = c("GP1", "GP2"),
    sample = c("S1", "S2"),
    value = c(10, 20)
  )

  res <- glyread:::.aggregate_long(tidy_df)

  expect_equal(nrow(res), 2L)
  expect_false("gp_id" %in% colnames(res))
  expect_equal(res$value, c(10, 20))
})
