#' Read Byonic-pGlycoQuant result
#'
#' If you used Byonic for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads in a pGlycoQuant result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @section Which file to use?:
#' You should use the "Quant.spectra.list" file in the pGlycoQuant result folder.
#' Files from Byonic result folder are not needed.
#' For instructions on how to use Byonic and pGlycoQuant, please refer to
#' the manual: [pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).
#'
#' @section Multisite glycopeptides:
#' Some glycopeptides can have more than one glycosylation site.
#' By default, they are expanded into multiple rows with the same quantification value
#' but different `protein_site` and `glycan_composition`.
#' Set `multisite = "drop"` to remove multisite glycopeptides instead.
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
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures (if `parse_structure = TRUE`).
#'
#' @inheritParams read_pglyco3_pglycoquant
#' @inheritParams read_byonic_byologic
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()]
#' @export
read_byonic_pglycoquant <- function(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db",
  parse_structure = TRUE,
  multisite = "expand"
) {
  # ----- Check arguments -----
  .validate_read_args(
    fp = fp,
    file_extensions = ".list",
    sample_info = sample_info,
    quant_method = quant_method,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    orgdb = orgdb
  )
  checkmate::assert_choice(multisite, c("expand", "drop"))

  # ----- Read data -----
  # Keep all PSM-level columns for variable info
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    df <- .read_byonic_df(fp)
    tidy_df <- .tidy_byonic_pglycoquant(df, orgdb, multisite)
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

.tidy_byonic_pglycoquant <- function(df, orgdb, multisite = "expand") {
  df %>%
    .handle_byonic_pglycoquant_multisite_rows(multisite) %>%
    .standardize_byonic_pglycoquant_columns() %>%
    .add_gene_symbols(orgdb) %>%
    .pivot_longer_pglycoquant()
}

.handle_byonic_pglycoquant_multisite_rows <- function(df, multisite) {
  switch(
    multisite,
    expand = .expand_byonic_pglycoquant_multisite_rows(df),
    drop = df %>%
      .drop_byonic_multisite_rows(
        glycan_col = "Composition",
        source = "Byonic-pGlycoQuant"
      ) %>%
      .expand_byonic_pglycoquant_multisite_rows()
  )
}

.read_byonic_df <- function(fp) {
  col_types <- readr::cols(
    Peptide = readr::col_character(),
    `Protein Name` = readr::col_character(),
    Position = readr::col_integer(),
    Composition = readr::col_character(),
  )
  suppressWarnings(
    suppressMessages(readr::read_tsv(
      fp,
      col_types = col_types,
      progress = FALSE
    )),
    classes = "vroom_mismatched_column_name"
  )
}

#' Expand Byonic-pGlycoQuant multisite glycopeptides
#'
#' pGlycoQuant reports Byonic multisite glycopeptides as one row with
#' comma-separated glycan compositions and multiple modified N residues in
#' `Peptide`. This helper expands those rows to one row per glycosylation site.
#'
#' @param df A Byonic-pGlycoQuant tibble.
#'
#' @returns A tibble with one row per glycosylation site.
#' @noRd
.expand_byonic_pglycoquant_multisite_rows <- function(df) {
  # Identify multisite glycopeptides (those with commas in Composition column)
  is_multisite <- stringr::str_detect(df$Composition, stringr::fixed(","))
  n_multisite <- sum(is_multisite)

  if (n_multisite > 0) {
    perc_multisite <- round(n_multisite / nrow(df) * 100, 1)
    cli::cli_alert_info(
      "Found {.val {n_multisite}} of {.val {nrow(df)}} ({.val {perc_multisite}}%) multisite PSMs. Expanding these entries to site-specific rows."
    )
  }

  df %>%
    dplyr::mutate(
      glycan_composition = purrr::map(
        .data$Composition,
        .split_byonic_glycan_compositions
      ),
      peptide_site = purrr::map(
        .data$Peptide,
        .extract_byonic_pglycoquant_glycosites
      ),
      n_glycans = purrr::map_int(.data$glycan_composition, length),
      n_sites = purrr::map_int(.data$peptide_site, length)
    ) %>%
    .check_byonic_glycan_site_pairing(
      source = "Byonic-pGlycoQuant",
      site_column = "Peptide"
    ) %>%
    dplyr::mutate(first_peptide_site = purrr::map_int(.data$peptide_site, min)) %>%
    tidyr::unnest_longer(c("glycan_composition", "peptide_site")) %>%
    dplyr::mutate(
      peptide_site = as.integer(.data$peptide_site),
      protein_site = as.integer(
        .data$Position + .data$peptide_site - .data$first_peptide_site
      )
    ) %>%
    dplyr::select(-all_of(c("first_peptide_site", "n_glycans", "n_sites")))
}

#' Standardize expanded Byonic-pGlycoQuant columns
#'
#' @param df A Byonic-pGlycoQuant tibble after multisite row expansion.
#'
#' @returns A wide-format tibble with standard glycopeptide variable columns
#'   and pGlycoQuant intensity columns.
#' @noRd
.standardize_byonic_pglycoquant_columns <- function(df) {
  var_cols <- c(
    "peptide",
    "peptide_site",
    "protein",
    "protein_site",
    "glycan_composition"
  )

  df %>%
    dplyr::mutate(
      peptide = .clean_byonic_pglycoquant_peptide(.data$Peptide),
      # >sp|P19652|A1AG2_HUMAN -> P19652, >sp|P08185-1|CBG_HUMAN -> P08185-1
      protein = .extract_uniprot_accession(.data$`Protein Name`)
    ) %>%
    dplyr::select(all_of(var_cols), tidyselect::starts_with("Intensity"))
}

#' Extract glycosylated peptide sites from a Byonic peptide
#'
#' @param peptide A Byonic peptide string with flanking residues.
#'
#' @returns An integer vector of peptide positions carrying glycans.
#' @noRd
.extract_byonic_pglycoquant_glycosites <- function(peptide) {
  if (is.na(peptide)) {
    return(NA_integer_)
  }

  peptide_core <- .extract_byonic_pglycoquant_peptide_core(peptide)
  glycan_matches <- stringr::str_locate_all(peptide_core, "N\\[[^\\]]+\\]")[[1]]

  if (nrow(glycan_matches) > 0) {
    return(purrr::map_int(
      glycan_matches[, "start"],
      \(match_start) {
        prefix <- stringr::str_sub(peptide_core, 1L, match_start)
        nchar(stringr::str_remove_all(prefix, "\\[[^\\]]+\\]"))
      }
    ))
  }

  lowercase_n_matches <- stringr::str_locate_all(peptide_core, "n")[[1]]
  as.integer(lowercase_n_matches[, "start"])
}

#' Remove flanking residues and modifications from a Byonic peptide
#'
#' @param peptide A Byonic peptide string with flanking residues.
#'
#' @returns A clean uppercase peptide sequence.
#' @noRd
.clean_byonic_pglycoquant_peptide <- function(peptide) {
  peptide %>%
    .extract_byonic_pglycoquant_peptide_core() %>%
    stringr::str_remove_all("\\[[^\\]]+\\]") %>%
    stringr::str_to_upper()
}

#' Extract the peptide core from a Byonic peptide
#'
#' @param peptide A Byonic peptide string with flanking residues.
#'
#' @returns The peptide sequence between flanking residues.
#' @noRd
.extract_byonic_pglycoquant_peptide_core <- function(peptide) {
  stringr::str_sub(peptide, 3L, -3L)
}
