test_that("read_glyhunter works", {
  exp <- suppressMessages(read_glyhunter(
    test_path("data/glyhunter-result.csv"),
    glycan_type = "N"
  ))

  expect_s4_class(exp, "GlycomicSE")
  expect_equal(.test_metadata(exp)$exp_type, "glycomics")
  expect_s3_class(.test_var_info(exp)$glycan_composition, "glyrepr_composition")
})

test_that("read_glyhunter works for NP preset", {
  exp <- suppressMessages(read_glyhunter(
    test_path("data/glyhunter-np-result.csv"),
    glycan_type = "N",
    preset = "Fu_NP_2026"
  ))

  expect_s4_class(exp, "GlycomicSE")
  expect_equal(.test_metadata(exp)$exp_type, "glycomics")
  expect_s3_class(.test_var_info(exp)$glycan_composition, "glyrepr_composition")
  expect_true(all(c("nL", "nE") %in% colnames(.test_var_info(exp))))
})

test_that("read_glyhunter works with sample name converter", {
  exp <- suppressMessages(read_glyhunter(
    test_path("data/glyhunter-result.csv"),
    glycan_type = "N",
    sample_name_converter = function(samples) {
      paste0("Sample_", samples)
    }
  ))

  expect_true(all(grepl("^Sample_", colnames(.test_expr_mat(exp)))))
})
