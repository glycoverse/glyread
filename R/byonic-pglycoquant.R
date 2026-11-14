#' Read Byonic-pGlycoQuant result
#'
#' If you used Byonic for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads in a pGlycoQuant result file and returns a [glyexp::GlycoproteomicSE()] object.
#' Currently only label-free quantification is supported.
#'
#' @section Which file to use?:
#' You should use the "Quant.spectra.list" file in the pGlycoQuant result folder.
#' Files from Byonic result folder are not needed.
#' For instructions on how to use Byonic and pGlycoQuant, please refer to
#' the manual: [pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).
#'
#' @section Glycosylation sites:
#' For each PSM, the reader matches the monoisotopic residue mass of each
#' reported glycan to a bracketed peptide modification within 0.001 Da.
#' Each peptide position can carry one reported glycan. All amino acids are
#' eligible. Protein positions are `Position + peptide_site - 1`, where
#' `Position` is the start of the peptide in the protein.
#'
#' The reader retains assignments shared by all maximum one-to-one matches.
#' Unresolved sites or site-specific compositions are recorded as `NA`,
#' with one warning summarizing the affected PSMs. If several unresolved
#' assignments give the same output feature, its PSM intensity is counted once.
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
#' @returns An [glyexp::GlycoproteomicSE()] object.
#' @seealso [glyexp::GlycoproteomicSE()], [glyrepr::glycan_composition()]
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
  fp <- .read_file_input(fp, ".list")
  .validate_read_args(
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

# Expand PSMs before sample intensities are pivoted or aggregated.
.expand_byonic_pglycoquant_multisite_rows <- function(df) {
  compositions <- lapply(df$Composition, .split_byonic_glycan_compositions)
  n_multisite <- sum(lengths(compositions) > 1L)
  if (n_multisite > 0L) {
    cli::cli_alert_info(
      "Expanding {n_multisite} multisite Byonic-pGlycoQuant PSM{?s} to site-specific rows."
    )
  }

  # Parse each distinct source composition once across all spectra.
  lookup <- .byonic_glycan_lookup(unique(as.character(unlist(compositions))))
  pairs <- Map(
    function(peptide, glycans) {
      ix <- match(glycans, lookup$source)
      .pair_byonic_glycans(peptide, lookup$composition[ix], lookup$mass[ix])
    },
    df$Peptide,
    compositions
  )

  unresolved <- sum(vapply(pairs, anyNA, logical(1)))
  if (unresolved > 0L) {
    cli::cli_warn(
      "Unresolved glycan-site assignments in {unresolved} of {nrow(df)} Byonic-pGlycoQuant PSMs. Unknown sites or compositions were retained as NA.",
      class = "glyread_unresolved_byonic_sites"
    )
  }

  df |>
    dplyr::mutate(pairs = unname(pairs)) |>
    tidyr::unnest("pairs", ptype = list(pairs = .byonic_site_pairs())) |>
    dplyr::mutate(protein_site = .data$Position + .data$peptide_site - 1L)
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
