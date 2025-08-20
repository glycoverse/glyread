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
#' Multisite glycopeptides are supported but their `protein_site` will be set to `NA` 
#' since the exact site of glycosylation cannot be determined unambiguously.
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
    .handle_multisite_byonic() %>%
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
    suppressMessages(readr::read_tsv(fp, col_types = col_types, progress = FALSE)),
    classes = "vroom_mismatched_column_name"
  )
}

.handle_multisite_byonic <- function(df) {
  # Identify multisite glycopeptides (those with commas in Composition column)
  is_multisite <- stringr::str_detect(df$Composition, stringr::fixed(","))
  n_multisite <- sum(is_multisite)
  
  if (n_multisite > 0) {
    perc_multisite <- round(n_multisite / nrow(df) * 100, 1)
    cli::cli_alert_info("Found {.val {n_multisite}} of {.val {nrow(df)}} ({.val {perc_multisite}}%) multisite PSMs. Setting protein_site to NA for these entries.")
  }
  
  # Keep all rows but mark multisite ones for special handling
  df$is_multisite <- is_multisite
  df
}

.refine_byonic_pglycoquant_columns <- function(df) {
  var_cols <- c("peptide", "peptide_site", "protein", "protein_site", "glycan_composition")
  df %>%
    .convert_byonic_columns() %>%
    dplyr::select(all_of(var_cols), tidyselect::starts_with("Intensity"))
}

.convert_byonic_columns <- function(df) {
  result <- df %>%
    dplyr::rename(
      peptide = "Peptide",
      protein = "Protein Name",
      protein_site = "Position",
      glycan_composition = "Composition"
    ) %>%
    dplyr::mutate(
      # K.N[+203.07937]GTR.G -> K.nGTR.G (mark glycosylated N)
      peptide = stringr::str_replace(.data$peptide, "N\\[.+?\\]", "n"),
      # K.nGTR.G -> nGTR (remove prefix and suffix)
      peptide = stringr::str_sub(.data$peptide, 3L, -3L),
      # Remove all other modifications like C[+57.02146] -> C
      peptide = stringr::str_remove_all(.data$peptide, "\\[.+?\\]"),
      # Extract peptide site (position of glycosylated residue)
      peptide_site = stringr::str_locate(.data$peptide, "n")[, "start"],
      # nGTR -> NGTR (restore N)
      peptide = stringr::str_replace(.data$peptide, "n", "N"),
      # >sp|P19652|A1AG2_HUMAN -> P19652, >sp|P08185-1|CBG_HUMAN -> P08185-1
      protein = .extract_uniprot_accession(.data$protein),
      # HexNAc(4)Hex(4)Fuc(1)NeuAc(1) -> HexNAc(4)Hex(4)dHex(1)NeuAc(1)
      glycan_composition = stringr::str_replace(.data$glycan_composition, "Fuc", "dHex")
    )
  
  # Set protein_site to NA for multisite glycopeptides (if is_multisite column exists)
  if ("is_multisite" %in% colnames(result)) {
    result <- dplyr::mutate(result, 
      protein_site = dplyr::if_else(.data$is_multisite, NA_integer_, .data$protein_site)
    )
  }
  
  result
}