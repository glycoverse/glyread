test_that("read_glyhunter works", {
  exp <- suppressMessages(read_glyhunter(
    test_path("glyhunter-result.csv"),
    glycan_type = "N"
  ))

  expect_s3_class(exp, "glyexp_experiment")
  expect_equal(exp$meta_data$exp_type, "glycomics")
  expect_s3_class(exp$var_info$glycan_composition, "glyrepr_composition")
})
