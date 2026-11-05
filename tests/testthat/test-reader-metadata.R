test_that("reader metadata checks allow additional constructor fields", {
  expected <- list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  )
  for (extra in list(
    list(),
    list(scale = "raw"),
    list(scale = "raw", future = TRUE)
  )) {
    se <- SummarizedExperiment::SummarizedExperiment(
      metadata = c(extra, expected)
    )
    .expect_reader_metadata(se, expected)
  }
})

test_that("reader metadata checks reject missing or incorrect required fields", {
  expected <- list(
    exp_type = "glycoproteomics",
    glycan_type = "N",
    quant_method = "label-free"
  )
  for (field in names(expected)) {
    missing <- c(expected, list(scale = "raw"))
    missing[[field]] <- NULL
    incorrect <- c(expected, list(scale = "raw"))
    incorrect[[field]] <- "incorrect"
    for (metadata in list(missing, incorrect)) {
      se <- SummarizedExperiment::SummarizedExperiment(metadata = metadata)
      expect_failure(.expect_reader_metadata(se, expected))
    }
  }
})
