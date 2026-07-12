test_that("numeric and factor sample identifiers are normalized", {
  numeric_info <- tibble::tibble(sample = c(1, 2))
  factor_info <- tibble::tibble(sample = factor(c("1", "2")))

  default_result <- suppressMessages(
    glyread:::.process_sample_info(NULL, c(1, 2))
  )
  numeric_result <- glyread:::.process_sample_info(
    numeric_info,
    c("1", "2")
  )
  factor_result <- glyread:::.process_sample_info(
    factor_info,
    c("1", "2")
  )

  expect_identical(default_result$sample, c("1", "2"))
  expect_identical(numeric_result$sample, c("1", "2"))
  expect_identical(factor_result$sample, c("1", "2"))
})
