#' Read Byonic-Byologic result
#'
#' If you used Byonic for intact glycopeptide identification,
#' and used Byologic for quantification, this is the function for you.
#' It reads in a result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @section Which file to use?:
#' Open the .blgc file in the result folder with PMI-Byos.
#' In the "Peptide List" panel (usually on the bottom right),
#' click "Export content of the table to a CSV file" button.
#' The exported .csv file is the file you should use.
#'
#' @section Multisite glycopeptides:
#' Currently, only single-site glycopeptides are supported.
#' Multisite glycopeptides will be removed.
#'
#' @inheritSection read_pglyco3_pglycoquant Sample information
#' @inheritSection read_pglyco3_pglycoquant Aggregation
#'
#' @section Variable information:
#' The following columns could be found in the variable information tibble:
#' - `peptide`: character, peptide sequence
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accession
#' - `protein_site`: integer, site of glycosylation on protein
#' - `gene`: character, gene name (symbol)
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#'
#' @inheritParams read_pglyco3_pglycoquant
#' @param orgdb name of the OrgDb package to use for UniProt to gene symbol conversion.
#'  Default is "org.Hs.eg.db".
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()]
#' @export
read_byonic_byologic <- function(
  fp,
  sample_info = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db"
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".csv")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)

  # ----- Read data -----
  # Keep all PSM-level columns for variable info
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    df <- .read_byonic_byologic_df(fp)
    tidy_df <- .tidy_byonic_byologic(df, orgdb)
    exp <- .read_template(
      tidy_df,
      sample_info,
      glycan_type,
      quant_method,
      sample_name_converter,
      composition_parser = glyrepr::as_glycan_composition
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.read_byonic_byologic_df <- function(fp) {
  df_raw <- readr::read_csv(fp, progress = FALSE, show_col_types = FALSE)
  df_clean <- df_raw |> janitor::clean_names()

  # Define expected column types based on cleaned names
  expected_cols <- readr::cols(
    row_number = readr::col_character(),      # Row#
    protein_name = readr::col_character(),    # Protein\nname or Protein\r\nname
    sequence = readr::col_character(),        # Sequence
    glycans = readr::col_character(),         # Glycans
    xic_area_summed = readr::col_double(),    # XIC area\nsummed or XIC area\r\nsummed
    ms_alias_name = readr::col_character(),   # MS\nAlias name or MS\r\nAlias name
    mod_summary = readr::col_character(),     # Mod.\nSummary or Mod.\r\nSummary
    start_aa = readr::col_integer()           # Start\nAA or Start\r\nAA
  )

  # Apply type conversion to the cleaned data
  df_typed <- df_clean |>
    dplyr::select(dplyr::any_of(names(expected_cols$cols))) |>
    readr::type_convert(col_types = expected_cols)

  return(df_typed)
}

.tidy_byonic_byologic <- function(df, orgdb) {
  df %>%
    .filter_byonic_byologic_rows() %>%
    .refine_byonic_byologic_columns() %>%
    .add_gene_symbols(orgdb)
}

.filter_byonic_byologic_rows <- function(df) {
  df %>%
    dplyr::filter(!is.na(.data$glycans)) %>%
    .collapse_byologic_rows() %>%
    .remove_multisite_byologic()
}

.remove_multisite_byologic <- function(df) {
  new_df <- dplyr::filter(df, !stringr::str_detect(.data$glycans, stringr::fixed(",")))
  n_removed <- nrow(df) - nrow(new_df)
  perc_removed <- round(n_removed / nrow(df) * 100, 1)
  cli::cli_alert_info("Removed {.val {n_removed}} of {.val {nrow(df)}} ({.val {perc_removed}}%) multisite glycopeptides.")
  new_df
}

# Collapse hierarchical search-result rows to peptide–sample abundances
#
# Many MS search engines export peptide tables where the *Row* column
# encodes a three-level hierarchy:
# * **Level 1** (“1”, “2”, …) — peptide aggregated across all samples.
# * **Level 2** (“1.1”, “1.2”, …) — peptide intensity in a single sample.
# * **Level 3** (“1.1.1”, “1.1.2”, …) — individual PSMs within that sample.
#
# This helper keeps exactly one row per **peptide × sample** pair by
# applying the following rules:
# 1. If a peptide appears in *multiple* samples, discard its level-1 row
#    and keep all level-2 rows.
# 2. If a peptide appears in *one* sample only, keep the single level-1
#    row (its intensity already equals that level-2 value).
# 3. Always drop level-3 rows (individual PSMs).
#
# The algorithm infers the hierarchy from the number of dots (`.`) in
# *Row*, so it is robust to any absolute numbering. The result can be
# returned in long format or pivoted to wide format for downstream
# statistics or visualization.
.collapse_byologic_rows <- function(df) {
  df %>%
    dplyr::mutate(
      depth = stringr::str_count(.data$row_number, stringr::fixed(".")) + 1,
      lvl1 = stringr::str_split_i(.data$row_number, stringr::fixed("."), 1L)
    ) %>%
    dplyr::mutate(has_lvl2 = any(.data$depth == 2), .by = "lvl1") %>%
    dplyr::filter(dplyr::if_else(.data$has_lvl2, .data$depth == 2, .data$depth == 1)) %>%
    dplyr::select(-all_of(c("depth", "lvl1", "has_lvl2")))
}

.refine_byonic_byologic_columns <- function(df) {
  selected_cols <- c("peptide", "peptide_site", "protein", "protein_site", "glycan_composition", "sample", "value")
  df %>%
    dplyr::mutate(
      # sp|P02765|FETUA_HUMAN... -> P02765
      protein = stringr::str_split_i(.data$protein_name, stringr::fixed("|"), 2L),
      # K.EHEGAIYPDnTTDFQR.A -> EHEGAIYPDnTTDFQR
      peptide = stringr::str_split_i(.data$sequence, stringr::fixed("."), 2L),
      # EHEGAIYPDnTTDFQR -> EHEGAIYPDNTTDFQR (n -> N)
      peptide = stringr::str_to_upper(.data$peptide),
      # Fuc -> dHex
      glycan_composition = stringr::str_replace(.data$glycans, "Fuc", "dHex"),
      # Add peptide_site
      peptide_site = stringr::str_extract(.data$mod_summary, "N(\\d+)\\(NGlycan", group = 1),
      peptide_site = as.integer(.data$peptide_site),
      # Add protein_site
      protein_site = as.integer(.data$start_aa + .data$peptide_site - 1L)
    ) %>%
    dplyr::rename(all_of(c(sample = "ms_alias_name", value = "xic_area_summed"))) %>%
    dplyr::select(all_of(selected_cols))
}
