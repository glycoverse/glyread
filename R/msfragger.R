#' Read MSFragger-Glyco result
#'
#' MSFragger-Glyco is a software for glycopeptide identification and quantification.
#' This function reads in the result file and returns a [glyexp::experiment()] object.
#'
#' This function uses the "psm.tsv" file in each sample folder.
#' Sample names are extracted from the file paths.
#' They are the parent of each "psm.tsv" file.
#' For example, "msfragger_result/H1/psm.tsv" will be named "H1".
#'
#' @param dp The directory path of the MSFragger-Glyco result folder.
#' @inheritParams read_pglyco3_pglycoquant
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()]
#' @export
read_msfragger <- function(
  dp,
  sample_info = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
) {
  # ----- Check arguments -----
  checkmate::assert_directory_exists(dp, access = "r")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)

  # ----- Read data -----
  if (quant_method == "label-free") {
    fps <- fs::dir_ls(dp, glob = "*/psm.tsv", recurse = TRUE)
    if (length(fps) == 0) {
      rlang::abort("No psm.tsv files found in the specified directory.")
    }
    exp <- .read_msfragger_label_free(fps, sample_info, glycan_type, sample_name_converter)
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.read_msfragger_label_free <- function(fps, sample_info, glycan_type, sample_name_converter) {
  # ----- Read data -----
  cli::cli_progress_step("Reading data")
  df <- .read_msfragger_file(fps, sample_name_converter) %>%
    .filter_unsure_msfragger_sites() %>%
    .convert_msfragger_columns()

  # ---- Aggregate PSMs to glycopeptides -----
  cli::cli_progress_step("Aggregating PSMs to glycopeptides")
  aggr_res <- .aggregate_long(df)
  var_info <- aggr_res$var_info
  expr_mat <- aggr_res$expr_mat
  rownames(expr_mat) <- var_info$variable
  sample_info <- .process_sample_info(sample_info, colnames(expr_mat), glycan_type)

  # ----- Parse glycan compositions -----
  cli::cli_progress_step("Parsing glycan compositions")
  var_info <- dplyr::mutate(var_info,
    glycan_composition = .parse_msfragger_glycan_composition(.data$glycan_composition),
  )

  # ----- Pack Experiment -----
  cli::cli_progress_step("Packing experiment")
  glyexp::experiment(
    expr_mat, sample_info, var_info,
    exp_type = "glycoproteomics",
    glycan_type = glycan_type,
    quant_method = "label-free"
  )
}

.read_msfragger_file <- function(fps, sample_name_converter) {
  col_types <- readr::cols(
    `Spectrum File` = readr::col_character(),
    `Peptide` = readr::col_character(),
    `Charge` = readr::col_integer(),
    `Protein Start` = readr::col_integer(),
    `Intensity` = readr::col_double(),
    `Total Glycan Composition` = readr::col_character(),
    `Best Positions` = readr::col_character(),
    `Number Best Positions` = readr::col_integer(),
    `Protein ID` = readr::col_character(),
    `Gene` = readr::col_character()
  )

  # Read files and ensure file column is always present
  suppressWarnings(
    suppressMessages(df <- readr::read_tsv(fps, col_types = col_types, progress = FALSE, id = "file")),
    classes = "vroom_mismatched_column_name"
  )

  df <- dplyr::mutate(df, sample = basename(dirname(.data$file)))
  if (!is.null(sample_name_converter)) {
    df <- dplyr::mutate(df, sample = sample_name_converter(.data$sample))
  }
  df
}

# Some PSM does not have definite glycosite assignment,
# and some have multiple glycosties.
# We filter them out.
.filter_unsure_msfragger_sites <- function(df) {
  new_df <- dplyr::filter(df, .data$`Number Best Positions` == 1)
  n_removed <- nrow(df) - nrow(new_df)
  perc_removed <- round(n_removed / nrow(df) * 100, 1)
  cli::cli_alert_info("Removed {.val {n_removed}} ({.val {perc_removed}}%) uncertain or multisite PSMs.")
  new_df
}

.convert_msfragger_columns <- function(df) {
  df %>%
    dplyr::select(
      peptide = "Peptide",
      charge = "Charge",
      protein_start = "Protein Start",
      glycan_composition = "Total Glycan Composition",
      best_positions = "Best Positions",
      protein = "Protein ID",
      gene = "Gene",
      sample = "sample",
      value = "Intensity"
    ) %>%
    dplyr::mutate(
      peptide_site = as.integer(stringr::str_sub(.data$best_positions, 2L, -1L)),
      protein_site = .data$protein_start + .data$peptide_site - 1L,
    )
}

.parse_msfragger_glycan_composition <- function(x) {
  x %>%
    stringr::str_split_i(" % ", 1L) %>%
    stringr::str_replace("Fuc", "dHex") %>%
    glyrepr::as_glycan_composition()
}
