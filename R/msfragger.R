#' Read MSFragger-Glyco result
#'
#' @description
#' MSFragger-Glyco is a software for glycopeptide identification and quantification.
#' This function reads in the result file and returns a [glyexp::experiment()] object.
#'
#' Because MSFragger performs quantification on the PSM level for each sample,
#' and outputs the result to separate files,
#' this function cannot read and merge multiple result files into a single [glyexp::experiment()] object.
#' See details to understand why.
#'
#' As a workaround, use this code snippet to read the result into a list of `glyexp::experiment()`` objects,
#' aggregate the experiments to the "glycoform" level using `glyclean::aggregate()`,
#' and then use `glyexp::merge()` (need glyexp >= 0.7.0) to merge them into one:
#'
#' ```r
#' library(tidyverse)
#' library(glyread)
#' library(glyclean)
#' library(glyexp)
#'
#' # Read the result into a list of experiments
#' files <- c("psm1.tsv", "psm2.tsv", "psm3.tsv")
#' exps <- map(files, read_msfragger)  # ignore sample_info for now
#'
#' # Aggregate the experiments to the "glycoform" level
#' exps_agg <- map(exps, aggregate, to_level = "gf")  # or other levels
#'
#' # Merge the experiments into one
#' exp_merged <- reduce(exps_agg, merge)
#'
#' # Add sample information
#' # Say you already have a sample information tibble (`sample_df`) with a "sample" column
#' exp_merged$sample_info <- exp_merged$sample_info |>
#'   left_join(sample_df, by = "sample")
#' ```
#'
#' Note that we override the sample information with the new one,
#' which is not a recommended practice.
#' However, the code here is harmless,
#' and you can safely use it until we find a better solution.
#'
#' @details
#' # Why not read multiple files?
#'
#' This is the first solution we thought of.
#' However, we found it is not compatible with the current design of `glycoverse` settings.
#' MSFragger-Glyco results are on the PSM level,
#' each raw sample resulting in a separate file.
#' To merge them into one expression matrix, with rows as PSMs and columns as samples,
#' we need to identify the PSM that appears in different samples.
#' This is not trivial, as an ion (with the same charge and modification state)
#' can generate multiple PSMs with different retention time and m/z.
#'
#' One possible solution is to aggregate each PSM result into a glycoform level.
#' That is, for a glycoform (unique combination of glycan, protein, and glycosite),
#' we sum up the quantifications of all PSMs that belong to this glycoform.
#' Then the data can be merged into one expression matrix,
#' with rows as glycoforms and columns as samples.
#' This is exactly what the about code snippet does.
#'
#' So why not just incorporate this approach into `glyread`?
#' Well, the job of `glyread` should be to read search results as they are.
#' It should not be responsible for any data preprocessing.
#' Aggregating data into "glycoform" level is a preprocessing step.
#' In fact, it is not even the first step
#' (aggregation after normalization and imputation is more preferred).
#' So we put this function in the `glyclean` package.
#'
#' Maybe in the future, MSFragger-Glyco has a merged output like pGlycoQuant.
#' Then everyone is happy.
#' But for now, we leave this to the user to do.
#' Hint: you can wrap the snippet above into a function yourself.
#'
#' @param fp File path of the MSFragger result file.
#' @inheritParams read_pglyco3_pglycoquant
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()]
#' @export
read_msfragger <- function(
  fp,
  sample_info = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".tsv")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)

  # ----- Read data -----
  if (quant_method == "label-free") {
    exp <- .read_msfragger_label_free(fp, sample_info, glycan_type, sample_name_converter)
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.read_msfragger_label_free <- function(fp, sample_info, glycan_type, sample_name_converter) {
  # ----- Read data -----
  cli::cli_progress_step("Reading data")
  df <- .read_msfragger_file(fp)
  df <- .filter_unsure_msfragger_sites(df)
  var_info <- .extract_msfragger_var_info(df)
  sample_name <- .get_msfragger_sample_name(df)
  if (!is.null(sample_name_converter)) {
    sample_name <- sample_name_converter(sample_name)
  }
  sample_info <- .process_sample_info(sample_info, sample_name, glycan_type)
  expr_mat <- .extract_msfragger_expr_mat(df)
  colnames(expr_mat) <- sample_name
  rownames(expr_mat) <- var_info$variable

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

.read_msfragger_file <- function(fp) {
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
  suppressWarnings(
    suppressMessages(readr::read_tsv(fp, col_types = col_types, progress = FALSE)),
    classes = "vroom_mismatched_column_name"
  )
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

.extract_msfragger_expr_mat <- function(df) {
  as.matrix(df$Intensity, ncol = 1)
}

.get_msfragger_sample_name <- function(df) {
  stringr::str_split_i(df$`Spectrum File`[[1]], "\\\\", -2L)
}

.extract_msfragger_var_info <- function(df) {
  new_names <- c(
    peptide = "Peptide",
    charge = "Charge",
    protein_start = "Protein Start",
    glycan_composition = "Total Glycan Composition",
    best_positions = "Best Positions",
    protein = "Protein ID",
    gene = "Gene"
  )
  df %>%
    dplyr::select(all_of(new_names)) %>%
    dplyr::mutate(
      peptide_site = as.integer(stringr::str_sub(.data$best_positions, 2L, -1L)),
      protein_site = .data$protein_start + .data$peptide_site - 1L,
    ) %>%
    dplyr::select(-all_of(c("protein_start", "best_positions"))) %>%
    dplyr::mutate(variable = stringr::str_c("PSM", dplyr::row_number()), .before = 1)
}

.parse_msfragger_glycan_composition <- function(x) {
  x %>%
    stringr::str_split_i(" % ", 1L) %>%
    stringr::str_replace("Fuc", "dHex") %>%
    glyrepr::as_glycan_composition()
}
