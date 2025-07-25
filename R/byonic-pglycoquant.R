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
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db",
  parse_structure = TRUE
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".list")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)

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
    .remove_multisite_byonic() %>%
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

.remove_multisite_byonic <- function(df) {
  new_df <- dplyr::filter(df, !stringr::str_detect(.data$Composition, stringr::fixed(",")))
  n_removed <- nrow(df) - nrow(new_df)
  perc_removed <- round(n_removed / nrow(df) * 100, 1)
  cli::cli_alert_info("Removed {.val {n_removed}} of {.val {nrow(df)}} ({.val {perc_removed}}%) multisite PSMs.")
  new_df
}

.refine_byonic_pglycoquant_columns <- function(df) {
  var_cols <- c("peptide", "peptide_site", "protein", "protein_site", "glycan_composition")
  df %>%
    .convert_byonic_columns() %>%
    dplyr::select(all_of(var_cols), tidyselect::starts_with("Intensity"))
}

.convert_byonic_columns <- function(df) {
  df %>%
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
      # >sp|P19652|A1AG2_HUMAN -> P19652
      protein = stringr::str_extract(.data$protein, "(?:sp|tr)\\|(\\w+)\\|.*", group = 1),
      # HexNAc(4)Hex(4)Fuc(1)NeuAc(1) -> HexNAc(4)Hex(4)dHex(1)NeuAc(1)
      glycan_composition = stringr::str_replace(.data$glycan_composition, "Fuc", "dHex"),
    )
}

.add_gene_symbols <- function(df, orgdb) {
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    cli::cli_alert_info("Package `{orgdb}` not installed. Skipping gene symbol conversion.")
    return(df)
  }

  # Get unique protein IDs (uniprot IDs)
  unique_proteins <- unique(df$protein)
  unique_proteins <- unique_proteins[!is.na(unique_proteins)]

  if (length(unique_proteins) == 0) {
    cli::cli_alert_info("No valid protein IDs found, skipping gene symbol conversion.")
    return(df)
  }

  # Try to convert uniprot IDs to gene symbols
  tryCatch(
    {
      org_db <- utils::getFromNamespace(orgdb, orgdb)

      # Use mapIds to convert UNIPROT to SYMBOL
      mapped_symbols <- suppressMessages(AnnotationDbi::mapIds(
        org_db,
        keys = unique_proteins,
        column = "SYMBOL",
        keytype = "UNIPROT",
        multiVals = "first"
      ))

      id_mapping <- tibble::tibble(
        UNIPROT = names(mapped_symbols),
        SYMBOL = mapped_symbols
      ) %>%
        tidyr::drop_na("SYMBOL")

      # Join with var_info to add gene symbols
      df <- dplyr::left_join(
        df,
        id_mapping,
        by = c("protein" = "UNIPROT")
      ) %>%
        dplyr::rename(gene = "SYMBOL")

      # Report conversion success
      n_converted <- sum(!is.na(unique(df$gene)))
      n_total <- length(unique(df$protein))
      perc_converted <- round(n_converted / n_total * 100, 1)
      cli::cli_alert_info(
        "Converted {.val {n_converted}} of {.val {n_total}} ({.val {perc_converted}}%) protein IDs to gene symbols."
      )
    },
    error = function(e) {
      cli::cli_alert_warning(
        "Failed to convert protein IDs to gene symbols: {.val {e$message}}"
      )
    }
  )

  df
}