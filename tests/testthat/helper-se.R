.test_expr_mat <- function(x) {
  SummarizedExperiment::assay(x, "abundance")
}

.test_var_info <- function(x) {
  tibble::as_tibble(
    SummarizedExperiment::rowData(x),
    rownames = "variable"
  )
}

.test_sample_info <- function(x) {
  tibble::as_tibble(
    SummarizedExperiment::colData(x),
    rownames = "sample"
  )
}

.test_metadata <- function(x) {
  S4Vectors::metadata(x)
}
