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

#' Drop multisite Byonic rows when requested
#'
#' @param df A Byonic result tibble.
#' @param glycan_col Name of the column containing glycan compositions.
#' @param source Human-readable source name for messages.
#'
#' @returns The input tibble with multisite rows removed.
#' @noRd
.drop_byonic_multisite_rows <- function(df, glycan_col, source) {
  is_multisite <- stringr::str_detect(df[[glycan_col]], stringr::fixed(","))
  n_multisite <- sum(is_multisite, na.rm = TRUE)

  if (n_multisite > 0) {
    perc_multisite <- round(n_multisite / nrow(df) * 100, 1)
    cli::cli_alert_info(
      "Dropping {.val {n_multisite}} of {.val {nrow(df)}} ({.val {perc_multisite}}%) multisite {source} rows."
    )
  }

  df %>%
    dplyr::mutate(.byonic_is_multisite = is_multisite) %>%
    dplyr::filter(!.data$.byonic_is_multisite) %>%
    dplyr::select(-all_of(".byonic_is_multisite"))
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
