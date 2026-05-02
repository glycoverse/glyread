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
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures (if `parse_structure = TRUE`).
#' - `gp_id`: character, glycopeptide ID before multisite row expansion.
#'
#' @inheritParams read_pglyco3_pglycoquant
#' @param orgdb name of the OrgDb package to use for UniProt to gene symbol conversion.
#'  Default is "org.Hs.eg.db".
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
  parse_structure = TRUE
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

  # ----- Read data -----
  # Keep all PSM-level columns for variable info
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    df <- .read_byonic_df(fp)
    tidy_df <- .tidy_byonic_pglycoquant(df, orgdb)
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

.tidy_byonic_pglycoquant <- function(df, orgdb) {
  df %>%
    .expand_byonic_pglycoquant_multisite_rows() %>%
    .refine_byonic_pglycoquant_columns() %>%
    .add_gene_symbols(orgdb) %>%
    .pivot_longer_pglycoquant()
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
#' `Peptide`. This helper expands those rows to one row per glycosylation site
#' while preserving a shared `gp_id`.
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
      gp_key = paste(
        .data$Peptide,
        .data$`Protein Name`,
        .data$Position,
        .data$Composition,
        sep = "\r"
      ),
      gp_id = paste0("BPGQ", dplyr::dense_rank(.data$gp_key)),
      glycan_composition = purrr::map(
        .data$Composition,
        .split_byonic_pglycoquant_glycans
      ),
      peptide_site = purrr::map(
        .data$Peptide,
        .extract_byonic_pglycoquant_glycosites
      ),
      first_peptide_site = purrr::map_int(.data$peptide_site, min),
      n_glycans = purrr::map_int(.data$glycan_composition, length),
      n_sites = purrr::map_int(.data$peptide_site, length)
    ) %>%
    .check_byonic_pglycoquant_glycosite_pairing() %>%
    tidyr::unnest_longer(c("glycan_composition", "peptide_site")) %>%
    dplyr::mutate(
      peptide_site = as.integer(.data$peptide_site),
      protein_site = as.integer(
        .data$Position + .data$peptide_site - .data$first_peptide_site
      )
    ) %>%
    dplyr::select(-all_of(c("gp_key", "first_peptide_site", "n_glycans", "n_sites")))
}

.refine_byonic_pglycoquant_columns <- function(df) {
  var_cols <- c(
    "gp_id",
    "peptide",
    "peptide_site",
    "protein",
    "protein_site",
    "glycan_composition"
  )
  df %>%
    .convert_byonic_columns() %>%
    dplyr::select(all_of(var_cols), tidyselect::starts_with("Intensity"))
}

.convert_byonic_columns <- function(df) {
  expanded_cols <- c("gp_id", "glycan_composition", "peptide_site", "protein_site")
  if (!all(expanded_cols %in% colnames(df))) {
    df <- .expand_byonic_pglycoquant_multisite_rows(df)
  }

  result <- df %>%
    dplyr::mutate(
      peptide = .clean_byonic_pglycoquant_peptide(.data$Peptide),
      # >sp|P19652|A1AG2_HUMAN -> P19652, >sp|P08185-1|CBG_HUMAN -> P08185-1
      protein = .extract_uniprot_accession(.data$`Protein Name`)
    )

  result
}

#' Split Byonic-pGlycoQuant glycan compositions by site
#'
#' @param glycans A single Byonic-pGlycoQuant glycan composition string.
#'
#' @returns A character vector with one glycan composition per site.
#' @noRd
.split_byonic_pglycoquant_glycans <- function(glycans) {
  if (is.na(glycans)) {
    return(NA_character_)
  }

  glycans %>%
    stringr::str_replace_all("Fuc", "dHex") %>%
    stringr::str_split(stringr::fixed(","), simplify = FALSE) %>%
    purrr::pluck(1) %>%
    stringr::str_trim()
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

#' Check glycan-site pairing in Byonic-pGlycoQuant rows
#'
#' @param df A Byonic-pGlycoQuant tibble with list columns for glycan
#'   compositions and peptide sites.
#'
#' @returns The input tibble, invisibly checked.
#' @noRd
.check_byonic_pglycoquant_glycosite_pairing <- function(df) {
  invalid_rows <- which(df$n_glycans != df$n_sites)
  if (length(invalid_rows) > 0) {
    rlang::abort(c(
      "Cannot pair Byonic-pGlycoQuant glycan compositions with glycosylation sites.",
      i = "The number of comma-separated glycans must match the number of glycosylated N sites in `Peptide`.",
      x = glue::glue("Problematic rows: {toString(invalid_rows)}")
    ))
  }

  df
}
