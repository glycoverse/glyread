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
#' Some glycopeptides can have more than one glycosylation site.
#' In this case, it will be expanded into multiple rows with the same quantification value
#' but different `protein_site` and `glycan_composition`.
#' This is generally fine if downstream analyses are done at the glycoform level.
#' A `gp_id` column is also added to uniquely identify each glycopeptide before row expansion,
#' so that the multiple rows can be mapped back to the original glycopeptide if needed.
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
#' - `gp_id`: character, glycopeptide ID before multisite row expansion.
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
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db"
) {
  # ----- Check arguments -----
  .validate_read_args(
    fp = fp,
    file_extensions = ".csv",
    sample_info = sample_info,
    quant_method = quant_method,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    orgdb = orgdb
  )

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
    row_number = readr::col_character(), # Row#
    protein_name = readr::col_character(), # Protein\nname or Protein\r\nname
    sequence = readr::col_character(), # Sequence
    glycans = readr::col_character(), # Glycans
    xic_area_summed = readr::col_double(), # XIC area\nsummed or XIC area\r\nsummed
    ms_alias_name = readr::col_character(), # MS\nAlias name or MS\r\nAlias name
    mod_summary = readr::col_character(), # Mod.\nSummary or Mod.\r\nSummary
    start_aa = readr::col_integer() # Start\nAA or Start\r\nAA
  )

  # Apply type conversion to the cleaned data
  df_clean |>
    dplyr::select(dplyr::any_of(names(expected_cols$cols))) |>
    readr::type_convert(col_types = expected_cols)
}

.tidy_byonic_byologic <- function(df, orgdb) {
  df %>%
    dplyr::filter(!is.na(.data$glycans)) %>%
    .collapse_byologic_rows() %>%
    .expand_byologic_multisite_rows() %>%
    .standardize_byologic_columns() %>%
    .add_gene_symbols(orgdb)
}

#' Expand Byologic multisite glycopeptides
#'
#' Byologic reports multisite glycopeptides as one row with comma-separated
#' glycan compositions and semicolon-separated NGlycan entries in `mod_summary`.
#' This helper expands those rows to one row per glycosylation site while
#' preserving a shared `gp_id`.
#'
#' @param df A collapsed Byologic tibble.
#'
#' @returns A tibble with one row per glycosylation site.
#' @noRd
.expand_byologic_multisite_rows <- function(df) {
  if (!"row_number" %in% colnames(df)) {
    df <- dplyr::mutate(df, row_number = as.character(dplyr::row_number()))
  }

  # Identify multisite glycopeptides (those with commas in glycans column)
  is_multisite <- stringr::str_detect(df$glycans, stringr::fixed(","))
  n_multisite <- sum(is_multisite)

  if (n_multisite > 0) {
    perc_multisite <- round(n_multisite / nrow(df) * 100, 1)
    cli::cli_alert_info(
      "Found {.val {n_multisite}} of {.val {nrow(df)}} ({.val {perc_multisite}}%) multisite glycopeptides. Expanding these entries to site-specific rows."
    )
  }

  df %>%
    dplyr::mutate(
      gp_source = stringr::str_split_i(.data$row_number, stringr::fixed("."), 1L),
      gp_id = paste0("BGP", dplyr::dense_rank(.data$gp_source)),
      glycan_composition = purrr::map(.data$glycans, .split_byologic_glycans),
      peptide_site = purrr::map(.data$mod_summary, .extract_byologic_glycosites),
      n_glycans = purrr::map_int(.data$glycan_composition, length),
      n_sites = purrr::map_int(.data$peptide_site, length)
    ) %>%
    .check_byologic_glycosite_pairing() %>%
    tidyr::unnest_longer(c("glycan_composition", "peptide_site")) %>%
    dplyr::mutate(peptide_site = as.integer(.data$peptide_site)) %>%
    dplyr::select(-all_of(c("gp_source", "n_glycans", "n_sites")))
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
    dplyr::filter(dplyr::if_else(
      .data$has_lvl2,
      .data$depth == 2,
      .data$depth == 1
    )) %>%
    dplyr::select(-all_of(c("depth", "lvl1", "has_lvl2")))
}

#' Standardize expanded Byologic columns
#'
#' @param df A Byologic tibble after multisite row expansion.
#'
#' @returns A long-format tibble with standard glycopeptide variable columns,
#'   `sample`, and `value`.
#' @noRd
.standardize_byologic_columns <- function(df) {
  selected_cols <- c(
    "gp_id",
    "peptide",
    "peptide_site",
    "protein",
    "protein_site",
    "glycan_composition",
    "sample",
    "value"
  )

  result <- df %>%
    dplyr::mutate(
      # sp|P02765|FETUA_HUMAN... -> P02765, sp|P08185-1|CBG_HUMAN -> P08185-1
      protein = .extract_uniprot_accession(.data$protein_name),
      # K.EHEGAIYPDnTTDFQR.A -> EHEGAIYPDnTTDFQR
      peptide = stringr::str_split_i(.data$sequence, stringr::fixed("."), 2L),
      # EHEGAIYPDnTTDFQR -> EHEGAIYPDNTTDFQR (n -> N)
      peptide = stringr::str_to_upper(.data$peptide),
      protein_site = as.integer(.data$start_aa + .data$peptide_site - 1L)
    )

  result %>%
    dplyr::rename(all_of(c(
      sample = "ms_alias_name",
      value = "xic_area_summed"
    ))) %>%
    dplyr::select(all_of(selected_cols))
}

#' Split Byologic glycan compositions by glycosylation site
#'
#' @param glycans A single Byologic glycan composition string.
#'
#' @returns A character vector with one glycan composition per site.
#' @noRd
.split_byologic_glycans <- function(glycans) {
  if (is.na(glycans)) {
    return(NA_character_)
  }

  glycans %>%
    stringr::str_replace_all("Fuc", "dHex") %>%
    stringr::str_split(stringr::fixed(","), simplify = FALSE) %>%
    purrr::pluck(1) %>%
    stringr::str_trim()
}

#' Extract Byologic glycosylation sites
#'
#' @param mod_summary A single Byologic modification summary string.
#'
#' @returns An integer vector with one peptide site per glycan.
#' @noRd
.extract_byologic_glycosites <- function(mod_summary) {
  if (is.na(mod_summary)) {
    return(NA_integer_)
  }

  matches <- stringr::str_match_all(
    mod_summary,
    "N(\\d+)\\(NGlycan"
  )[[1]]

  as.integer(matches[, 2])
}

#' Check glycan-site pairing in Byologic rows
#'
#' @param df A collapsed Byologic tibble with list columns for glycan
#'   compositions and peptide sites.
#'
#' @returns The input tibble, invisibly checked.
#' @noRd
.check_byologic_glycosite_pairing <- function(df) {
  invalid_rows <- df$row_number[df$n_glycans != df$n_sites]
  if (length(invalid_rows) > 0) {
    rlang::abort(c(
      "Cannot pair Byologic glycan compositions with glycosylation sites.",
      i = "The number of comma-separated glycans must match the number of NGlycan sites in `mod_summary`.",
      x = glue::glue("Problematic rows: {toString(invalid_rows)}")
    ))
  }

  df
}
