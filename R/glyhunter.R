#' Read GlyHunter result
#'
#' This function reads in a [GlyHunter](https://github.com/fubin1999/glyhunter)
#' result file and returns a [glyexp::experiment()] object.
#'
#' @details
#' # Which file to use?
#'
#' Use the "summary_area.csv" file in the GlyHunter result folder.
#' No edit is necessary.
#'
#' # Variable information
#'
#' The following columns could be found in the variable information tibble:
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#'
#' If "Fu_NP_2026" preset is used, two additional columns will be included:
#' - `nL`: Number of α2,3-linked sialic acids.
#' - `nE`: Number of α2,6-linked sialic acids.
#'
#' @inheritSection read_pglyco3_pglycoquant Sample information
#'
#' @inheritParams read_pglyco3_pglycoquant
#' @param fp File path of the GlyHunter result file.
#' @param preset "Fu_NC_2026" (default) or "Fu_NP_2026".
#'   Use "Fu_NC_2026" to use the configuration following "DOI: 10.1038/s41467-026-68579-x".
#'   Use "Fu_NP_2026" to use the configuration in an unpublished Nature Protocols paper.
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()]
#' @export
read_glyhunter <- function(
  fp,
  sample_info = NULL,
  preset = c("Fu_NC_2026", "Fu_NP_2026"),
  glycan_type = "N",
  sample_name_converter = NULL
) {
  # ----- Check arguments -----
  .validate_read_args(
    fp = fp,
    file_extensions = ".csv",
    sample_info = sample_info,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter
  )
  preset <- rlang::arg_match(preset)
  # ----- Read data -----
  switch(
    preset,
    "Fu_NC_2026" = .read_glyhunter_nc(
      fp,
      sample_info,
      glycan_type,
      sample_name_converter
    ),
    "Fu_NP_2026" = .read_glyhunter_np(
      fp,
      sample_info,
      glycan_type,
      sample_name_converter
    )
  )
}

.read_glyhunter_nc <- function(
  fp,
  sample_info = NULL,
  glycan_type = "N",
  sample_name_converter = NULL
) {
  df <- suppressMessages(readr::read_csv(fp))

  # Rename samples
  if (!is.null(sample_name_converter)) {
    samples <- setdiff(colnames(df), "glycan")
    new_samples <- sample_name_converter(samples)
    colnames(df) <- c("glycan", new_samples)
  }

  # Prepare sample info
  samples <- setdiff(colnames(df), "glycan")
  sample_info <- .process_sample_info(sample_info, samples, glycan_type)

  # Prepare variable info
  var_info <- df |>
    dplyr::select(all_of(c("variable" = "glycan"))) |>
    dplyr::mutate(
      variable = stringr::str_remove_all(.data$variable, "\\[.*?\\]"),
      glycan_composition = glyrepr::as_glycan_composition(.data$variable),
      variable = .data$variable |>
        stringr::str_replace(stringr::fixed("Hex("), "H") |>
        stringr::str_replace(stringr::fixed("HexNAc("), "N") |>
        stringr::str_replace(stringr::fixed("dHex("), "F") |>
        stringr::str_replace(stringr::fixed("NeuAc("), "S") |>
        stringr::str_remove_all(stringr::fixed(")")),
    )

  # Prepare expression matrix
  expr_mat <- df |>
    tibble::column_to_rownames("glycan") |>
    as.matrix()
  rownames(expr_mat) <- var_info$variable

  # Pack experiment
  glyexp::experiment(
    expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = glycan_type
  )
}

.read_glyhunter_np <- function(
  fp,
  sample_info = NULL,
  glycan_type = "N",
  sample_name_converter = NULL
) {
  df <- suppressMessages(readr::read_csv(fp))

  # Rename samples
  if (!is.null(sample_name_converter)) {
    samples <- setdiff(colnames(df), "glycan")
    new_samples <- sample_name_converter(samples)
    colnames(df) <- c("glycan", new_samples)
  }

  # Prepare sample info
  samples <- setdiff(colnames(df), "glycan")
  sample_info <- .process_sample_info(sample_info, samples, glycan_type)

  # Prepare variable info
  var_info <- df |>
    dplyr::select(all_of(c("variable" = "glycan"))) |>
    dplyr::mutate(
      variable = .data$variable |>
        stringr::str_replace(stringr::fixed("Hex("), "H") |>
        stringr::str_replace(stringr::fixed("HexNAc("), "N") |>
        stringr::str_replace(stringr::fixed("dHex("), "F") |>
        stringr::str_replace("NeuAc\\[\\+13.*?\\]\\(", "L") |>
        stringr::str_replace("NeuAc\\[\\+16.*?\\]\\(", "E") |>
        stringr::str_remove_all(stringr::fixed(")")),
      nL = stringr::str_extract(.data$variable, "L(\\d+)", group = 1),
      nL = dplyr::if_else(is.na(.data$nL), 0L, as.integer(.data$nL)),
      nE = stringr::str_extract(.data$variable, "E(\\d+)", group = 1),
      nE = dplyr::if_else(is.na(.data$nE), 0L, as.integer(.data$nE)),
      nS = .data$nL + .data$nE,
      glycan_composition = stringr::str_remove_all(.data$variable, "[EL]\\d+"),
      glycan_composition = dplyr::if_else(
        .data$nS > 0,
        paste0(.data$glycan_composition, "S", .data$nS),
        .data$glycan_composition
      ),
      glycan_composition = glyrepr::as_glycan_composition(
        .data$glycan_composition
      )
    ) |>
    dplyr::select(-all_of("nS"))

  # Prepare expression matrix
  expr_mat <- df |>
    tibble::column_to_rownames("glycan") |>
    as.matrix()
  rownames(expr_mat) <- var_info$variable

  # Pack experiment
  glyexp::experiment(
    expr_mat,
    sample_info = sample_info,
    var_info = var_info,
    exp_type = "glycomics",
    glycan_type = glycan_type
  )
}
