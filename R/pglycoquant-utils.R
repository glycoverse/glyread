.pivot_longer_pglycoquant <- function(df) {
  df %>%
    tidyr::pivot_longer(
      tidyselect::starts_with("Intensity"),
      names_to = "sample",
      values_to = "value",
      names_pattern = "Intensity\\((.*)\\)"
    )
}