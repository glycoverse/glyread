# The template reading function performing the following steps:
# 1. Aggregating PSMs to glycopeptides
# 2. Sample name conversion
# 3. Process sample information
# 4. Extracting variable information and expression matrix
# 5. Parse glycan compositions and structures
# 6. Packing an experiment
#
# `tidy_df` should be a long-format tibble with the following columns:
# "sample", "value", and all other columns being variable information.
.read_template <- function(
  tidy_df,
  sample_info_arg,
  glycan_type,
  quant_method,
  sample_name_converter = NULL,
  composition_parser = NULL,
  structure_parser = NULL
) {
  # ----- 1. Aggregate PSMs to glycopeptides -----
  tidy_df <- .aggregate_long(tidy_df)

  # ----- 2. Sample name conversion -----
  if (!is.null(sample_name_converter)) {
    tidy_df <- dplyr::mutate(tidy_df, sample = sample_name_converter(.data$sample))
  }

  # ----- 3. Process sample information -----
  samples <- unique(tidy_df$sample)
  sample_info <- .process_sample_info(sample_info_arg, samples, glycan_type)

  # ----- 4. Extracting variable information and expression matrix -----
  wide_df <- tidy_df %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "value")
  expr_mat <- as.matrix(wide_df[, samples])
  var_info <- wide_df %>%
    dplyr::select(-all_of(samples)) %>%
    dplyr::mutate(variable = paste0("GP", dplyr::row_number()), .before = 1)
  rownames(expr_mat) <- var_info$variable

  # ----- 5. Parse glycan compositions and structures -----
  cli::cli_progress_step("Parsing glycan compositions and structures")
  var_info <- dplyr::mutate(
    var_info,
    glycan_composition = composition_parser(.data$glycan_composition)
  )
  if (!is.null(structure_parser)) {
    var_info <- dplyr::mutate(
      var_info,
      glycan_structure = structure_parser(.data$glycan_structure)
    )
  }

  # ----- 6. Packing Experiment -----
  glyexp::experiment(
    expr_mat, sample_info, var_info,
    exp_type = "glycoproteomics",
    glycan_type = glycan_type,
    quant_method = quant_method
  )
}