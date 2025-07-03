# Extract expression matrix from pGlycoQuant result
# Note: rownames and colnames are not set here
.extract_expr_mat_from_pglycoquant <- function(df) {
  expr_mat <- df %>%
    dplyr::select(tidyselect::starts_with("Intensity")) %>%
    as.matrix()
  expr_mat[expr_mat == 0] <- NA
  colnames(expr_mat) <- NULL
  expr_mat
}


.extract_sample_names_from_pglycoquant <- function(df, sample_name_converter) {
  intensity_cols <- colnames(df)[stringr::str_starts(colnames(df), "Intensity")]
  samples <- stringr::str_extract(intensity_cols, "Intensity\\((.*)\\)", group = 1)
  if (!is.null(sample_name_converter)) {
    new_samples <- sample_name_converter(samples)
    if (length(new_samples) != length(samples)) {
      rlang::abort("Sample name converter must return the same number of samples.")
    }
    samples <- new_samples
  }
  samples
}
