# Aggregate PSMs to glycopeptides by summing expression values
# for unique combinations of glycopeptide-defining columns

# `df` should have the following columns:
# 1. the aggregation columns, 2. "sample", 3. "value"
.aggregate_long <- function(df) {
  aggr_cols <- c(
    "peptide",
    "peptide_site",
    "protein",
    "protein_site",
    "gene",
    "glycan_composition",
    "glycan_structure",
    "proteins",
    "genes",
    "protein_sites"
  )

  group_cols <- intersect(aggr_cols, colnames(df))
  value_df <- df %>%
    dplyr::summarise(
      value = sum(.data$value, na.rm = TRUE),
      .by = any_of(c(group_cols, "sample"))
    )

  if ("gp_id" %in% colnames(df)) {
    gp_id_df <- df %>%
      dplyr::summarise(
        gp_id = .collapse_trace_ids(.data$gp_id),
        .by = any_of(group_cols)
      )
    value_df <- dplyr::left_join(value_df, gp_id_df, by = group_cols)
  }

  value_df %>%
    dplyr::mutate(value = dplyr::if_else(.data$value == 0, NA, .data$value))
}

#' Collapse trace IDs while preserving first-seen order
#'
#' @param x A character vector of trace IDs.
#'
#' @returns A semicolon-separated character scalar, or `NA_character_` if all
#'   IDs are missing.
#' @noRd
.collapse_trace_ids <- function(x) {
  x <- unique(stats::na.omit(x))
  if (length(x) == 0) {
    return(NA_character_)
  }

  paste(x, collapse = ";")
}
