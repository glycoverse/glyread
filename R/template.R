# The template reading function performing the following steps:
# 1. Aggregating PSMs to glycopeptides
# 2. Sample name conversion
# 3. Process sample information
# 4. Extracting variable information and expression matrix
# 5. Parse glycan compositions and structures
# 6. Packing a GlycoproteomicSE
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
  structure_parser = NULL,
  parse_structure = TRUE
) {
  # ----- 1. Aggregate PSMs to glycopeptides -----
  tidy_df <- .aggregate_long(tidy_df)

  # ----- 2. Sample name conversion -----
  if (!is.null(sample_name_converter)) {
    tidy_df <- dplyr::mutate(
      tidy_df,
      sample = sample_name_converter(.data$sample)
    )
  }

  # ----- 3. Process sample information -----
  samples <- unique(tidy_df$sample)
  sample_info <- .process_sample_info(sample_info_arg, samples)

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
  if (parse_structure && !is.null(structure_parser)) {
    var_info <- dplyr::mutate(
      var_info,
      glycan_structure = structure_parser(.data$glycan_structure)
    )
  } else if (!parse_structure && "glycan_structure" %in% colnames(var_info)) {
    # Remove glycan_structure column if parse_structure is FALSE
    var_info <- dplyr::select(var_info, -any_of("glycan_structure"))
  }

  # ----- 6. Pack GlycoproteomicSE -----
  .new_glycoproteomic_se(
    expr_mat,
    sample_info,
    var_info,
    glycan_type = glycan_type,
    quant_method = quant_method
  )
}

#' Create a GlycoproteomicSE from reader output
#'
#' @param expr_mat A numeric abundance matrix.
#' @param sample_info A sample information data frame with a `sample` column.
#' @param var_info A variable information data frame with a `variable` column.
#' @param glycan_type The glycan type.
#' @param quant_method The quantification method.
#'
#' @returns A [glyexp::GlycoproteomicSE()] object.
#' @noRd
.new_glycoproteomic_se <- function(
  expr_mat,
  sample_info,
  var_info,
  glycan_type,
  quant_method
) {
  ids <- .prepare_se_ids(expr_mat, sample_info, var_info)
  expr_mat <- ids$expr_mat
  sample_info <- ids$sample_info
  var_info <- ids$var_info

  var_info$variable <- .standardize_glycoproteomic_variables(var_info)
  rownames(expr_mat) <- var_info$variable

  .muffle_all_missing_assay_warning(
    expr_mat,
    glyexp::GlycoproteomicSE(
      expr_mat,
      colData = .as_se_data(sample_info, "sample"),
      rowData = .as_se_data(var_info, "variable"),
      metadata = list(
        exp_type = "glycoproteomics",
        glycan_type = glycan_type,
        quant_method = quant_method
      )
    )
  )
}

#' Create a GlycomicSE from reader output
#'
#' @inheritParams .new_glycoproteomic_se
#'
#' @returns A [glyexp::GlycomicSE()] object.
#' @noRd
.new_glycomic_se <- function(expr_mat, sample_info, var_info, glycan_type) {
  ids <- .prepare_se_ids(expr_mat, sample_info, var_info)

  .muffle_all_missing_assay_warning(
    ids$expr_mat,
    glyexp::GlycomicSE(
      ids$expr_mat,
      colData = .as_se_data(ids$sample_info, "sample"),
      rowData = .as_se_data(ids$var_info, "variable"),
      metadata = list(
        exp_type = "glycomics",
        glycan_type = glycan_type
      )
    )
  )
}

#' Muffle the empty-min warning for an all-missing assay
#'
#' The glyco SE validity methods use `min(..., na.rm = TRUE)` to reject
#' negative abundance. Base R warns when the assay contains no non-missing
#' values, even though such an assay is valid.
#'
#' @param abundance A numeric abundance matrix.
#' @param code Code that constructs a glyco SE object.
#'
#' @returns The result of `code`.
#' @noRd
.muffle_all_missing_assay_warning <- function(abundance, code) {
  all_missing <- all(is.na(abundance))

  withCallingHandlers(
    code,
    warning = function(cnd) {
      is_empty_min_warning <- startsWith(
        conditionMessage(cnd),
        "no non-missing arguments to min"
      )
      if (all_missing && is_empty_min_warning) {
        rlang::cnd_muffle(cnd)
      }
    }
  )
}

#' Validate and align assay and dimension metadata identifiers
#'
#' @inheritParams .new_glycoproteomic_se
#'
#' @returns A list containing the aligned assay, sample information, and
#'   variable information.
#' @noRd
.prepare_se_ids <- function(expr_mat, sample_info, var_info) {
  checkmate::assert_names(colnames(sample_info), must.include = "sample")
  checkmate::assert_names(colnames(var_info), must.include = "variable")
  checkmate::assert_character(
    sample_info$sample,
    any.missing = FALSE,
    unique = TRUE
  )
  checkmate::assert_character(
    var_info$variable,
    any.missing = FALSE,
    unique = TRUE
  )

  if (!setequal(colnames(expr_mat), sample_info$sample)) {
    cli::cli_abort(
      "Sample identifiers in {.arg sample_info} must match assay column names."
    )
  }
  if (!setequal(rownames(expr_mat), var_info$variable)) {
    cli::cli_abort(
      "Variable identifiers in {.arg var_info} must match assay row names."
    )
  }

  list(
    expr_mat = expr_mat[var_info$variable, sample_info$sample, drop = FALSE],
    sample_info = sample_info,
    var_info = var_info
  )
}

#' Convert identifier-bearing metadata to an S4 DataFrame
#'
#' @param x A data frame.
#' @param id_col The identifier column to move to row names.
#'
#' @returns An [S4Vectors::DataFrame()] without the identifier column.
#' @noRd
.as_se_data <- function(x, id_col) {
  ids <- x[[id_col]]
  x[[id_col]] <- NULL
  S4Vectors::DataFrame(x, row.names = ids)
}

#' Standardize glycoproteomics variable identifiers
#'
#' @param var_info A variable information data frame.
#'
#' @returns A unique character vector of variable identifiers.
#' @noRd
.standardize_glycoproteomic_variables <- function(var_info) {
  variables <- paste(
    var_info$protein,
    dplyr::if_else(
      is.na(var_info$protein_site),
      "X",
      as.character(var_info$protein_site)
    ),
    as.character(var_info$glycan_composition),
    sep = "-"
  )

  tibble::tibble(variable = variables) |>
    dplyr::group_by(.data$variable) |>
    dplyr::mutate(
      variable = if (dplyr::n() > 1L) {
        paste0(.data$variable, "-", dplyr::row_number())
      } else {
        .data$variable
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::pull(.data$variable)
}
