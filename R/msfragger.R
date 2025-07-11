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
    cli::cli_progress_step("Reading data")
    fps <- fs::dir_ls(dp, glob = "*/psm.tsv", recurse = TRUE)
    if (length(fps) == 0) {
      rlang::abort("No psm.tsv files found in the specified directory.")
    }
    df <- .read_msfragger_file(fps)
    tidy_df <- .tidy_msfragger(df)
    exp <- .read_template(
      tidy_df,
      sample_info,
      glycan_type,
      quant_method,
      sample_name_converter,
      composition_parser = .parse_msfragger_glycan_composition
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.read_msfragger_file <- function(fps) {
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
  df
}

.tidy_msfragger <- function(df) {
  df %>%
    .filter_unsure_msfragger_sites() %>%
    .convert_msfragger_columns()
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
