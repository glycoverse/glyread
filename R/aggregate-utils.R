# Aggregate PSMs to glycopeptides by summing expression values
# for unique combinations of glycopeptide-defining columns
.aggregate_psms_to_glycopeptides <- function(var_info, expr_mat) {
  # Define glycopeptide-defining columns (after protein inference)
  aggr_cols <- c(
    "peptide", "peptide_site", "protein", "protein_site",
    "gene", "glycan_composition", "glycan_structure"
  )
  var_info <- dplyr::select(var_info, any_of(aggr_cols))
  sample_names <- colnames(expr_mat)
  res_df <- dplyr::bind_cols(var_info, tibble::as_tibble(expr_mat)) %>%
    tidyr::pivot_longer(all_of(sample_names), names_to = "sample", values_to = "value") %>%
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
