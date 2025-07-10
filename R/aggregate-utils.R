# Aggregate PSMs to glycopeptides by summing expression values
# for unique combinations of glycopeptide-defining columns

# `df` should have the following columns:
# 1. the aggregation columns, 2. "sample", 3. "value"
.aggregate_long <- function(df) {
  aggr_cols <- c(
    "peptide", "peptide_site", "protein", "protein_site",
    "gene", "glycan_composition", "glycan_structure",
    "proteins", "genes", "protein_sites"  # This three is for pGlyco3
  )
  sample_names <- unique(df$sample)
  res_df <- df %>%
    dplyr::summarise(value = sum(.data$value, na.rm = TRUE), .by = any_of(c(aggr_cols, "sample"))) %>%
    dplyr::mutate(value = dplyr::if_else(.data$value == 0, NA, .data$value)) %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "value")
  new_var_info <- res_df %>%
    dplyr::select(any_of(aggr_cols)) %>%
    dplyr::mutate(variable = paste0("GP", dplyr::row_number()), .before = 1)
  new_expr_mat <- as.matrix(res_df[, sample_names])
  rownames(new_expr_mat) <- new_var_info$variable
  list(var_info = new_var_info, expr_mat = new_expr_mat)
}

# This function is for pGlycoQuant results,
# in which extracting var_info and expr_mat directly from the raw result is more straighforward.
# So in those functions, like `read_pglyco3_pglycoquant()` and `read_byonic_pglycoquant()`,
# we first extract var_info and expr_mat, and then call this function to aggregate them.
.aggregate_wide <- function(var_info, expr_mat) {
  dplyr::bind_cols(var_info, tibble::as_tibble(expr_mat)) %>%
    tidyr::pivot_longer(all_of(colnames(expr_mat)), names_to = "sample", values_to = "value") %>%
    .aggregate_long()
}
