#' Split Byonic glycan compositions by glycosylation site
#'
#' @param glycans A single comma-separated Byonic glycan composition string.
#'
#' @returns A character vector with one glycan composition per site.
#' @noRd
.split_byonic_glycan_compositions <- function(glycans) {
  if (is.na(glycans)) {
    return(NA_character_)
  }

  glycans %>%
    stringr::str_replace_all("Fuc", "dHex") %>%
    stringr::str_split(stringr::fixed(","), simplify = FALSE) %>%
    purrr::pluck(1) %>%
    stringr::str_trim()
}

#' Check glycan-site pairing in expanded Byonic rows
#'
#' @param df A tibble with `n_glycans` and `n_sites` columns.
#' @param source Human-readable source name for error messages.
#' @param site_column Source column that carries the glycosylation sites.
#' @param row_id Row identifiers to report on mismatch.
#'
#' @returns The input tibble.
#' @noRd
.check_byonic_glycan_site_pairing <- function(
  df,
  source,
  site_column,
  row_id = seq_len(nrow(df))
) {
  invalid_rows <- row_id[df$n_glycans != df$n_sites]
  if (length(invalid_rows) > 0) {
    cli::cli_abort(c(
      "Cannot pair {source} glycan compositions with glycosylation sites.",
      i = "The number of comma-separated glycans must match the number of glycosylation sites in {.field {site_column}}.",
      x = "Problematic rows: {.val {toString(invalid_rows)}}"
    ))
  }

  df
}
