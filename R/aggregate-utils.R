# Aggregate PSMs to glycopeptides by summing expression values
# for unique combinations of glycopeptide-defining columns

# `df` should have the following columns:
# 1. the aggregation columns, 2. "sample", 3. "value"
.aggregate_long <- function(df) {
  aggr_cols <- c(
    "peptide", "peptide_site", "protein", "protein_site",
    "gene", "glycan_composition", "glycan_structure",
    "proteins", "genes", "protein_sites"
  )
  sample_names <- unique(df$sample)
  res_df <- df %>%
    dplyr::summarise(
      value = sum(.data$value, na.rm = TRUE), 
      .by = any_of(c(aggr_cols, "sample"))
    ) %>%
    dplyr::mutate(value = dplyr::if_else(.data$value == 0, NA, .data$value))
}
