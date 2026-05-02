test_that("gp_id is retained as trace metadata without splitting aggregation", {
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
  expect_equal(res$gp_id, "GP1;GP2")
})

test_that("gp_id trace metadata is stable across samples", {
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
  expect_equal(unique(res$gp_id), "GP1;GP2")
})
