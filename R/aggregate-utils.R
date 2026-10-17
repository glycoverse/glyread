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

  result <- df %>%
    dplyr::summarise(
      value = sum(.data$value, na.rm = TRUE),
      .by = any_of(c(aggr_cols, "sample"))
    ) %>%
    dplyr::mutate(value = dplyr::if_else(.data$value == 0, NA, .data$value))
  if ("protein_group" %in% names(df)) {
    keys <- intersect(aggr_cols, names(df))
    result <- dplyr::left_join(
      result,
      .protein_group_annotations(df, keys),
      by = keys
    )
  }
  result
}

# Combine source protein accessions without changing feature or sample grouping.
.protein_group_annotations <- function(df, keys) {
  df |>
    dplyr::summarise(
      protein_group = {
        groups <- .data$protein_group[!is.na(.data$protein_group)]
        proteins <- unique(unlist(strsplit(groups, ";", fixed = TRUE)))
        proteins <- proteins[nzchar(proteins)]
        if (length(proteins)) paste(proteins, collapse = ";") else NA_character_
      },
      .by = tidyselect::all_of(keys)
    )
}
