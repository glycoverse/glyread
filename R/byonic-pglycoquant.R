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
#'
#' @section Output:
#' This function returns a [glyexp::experiment()] object.
#'
#' The following columns could be found in the variable information tibble:
#' - `peptide`: character, peptide sequence
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accessions
#' - `protein_site`: integer, site of glycosylation on protein
#' - `gene`: character, gene symbols (if clusterProfiler package is available)
#'
#' @inheritParams read_pglyco3_pglycoquant
#' @param bitr_db name of the OrgDb package to use for UniProt to gene symbol conversion.
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
  bitr_db = "org.Hs.eg.db"
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
    exp <- .read_byonic_pglycoquant_label_free(
      fp, sample_info, glycan_type, sample_name_converter, bitr_db
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}


.read_byonic_pglycoquant_label_free <- function(
  fp,
  sample_info,
  glycan_type,
  sample_name_converter,
  bitr_db
) {
  # ----- Read data -----
  cli::cli_progress_step("Reading data")
  df <- .read_byonic_file_into_tibble(fp)
  df <- .remove_multisite_byonic(df)
  expr_mat <- .extract_expr_mat_from_pglycoquant(df)
  samples <- .extract_sample_names_from_pglycoquant(df, sample_name_converter)
  sample_info <- .process_sample_info(sample_info, samples, glycan_type)
  var_info <- .extract_var_info_from_byonic(df, bitr_db)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable

  # ----- Parse glycan compositions -----
  cli::cli_progress_step("Parsing glycan compositions")
  var_info <- dplyr::mutate(var_info,
    glycan_composition = glyrepr::as_glycan_composition(.data$glycan_composition)
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


.read_byonic_file_into_tibble <- function(fp) {
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


.extract_var_info_from_byonic <- function(df, bitr_db) {
  new_names <- c(
    peptide = "Peptide",
    protein = "Protein Name",
    protein_site = "Position",
    glycan_composition = "Composition"
  )

  var_info <- df %>%
    dplyr::select(all_of(new_names)) %>%
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
    ) %>%
    dplyr::mutate(variable = stringr::str_c("PSM", dplyr::row_number()), .before = 1) %>%
    .add_gene_symbols(bitr_db)

  var_info
}

.add_gene_symbols <- function(var_info, bitr_db) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE) || !requireNamespace(bitr_db, quietly = TRUE)) {
    cli::cli_alert_info("Package `clusterProfiler` and/or `{bitr_db}` not installed. Skipping gene symbol conversion.")
    return(var_info)
  }

  # Get unique protein IDs (uniprot IDs)
  unique_proteins <- unique(var_info$protein)
  unique_proteins <- unique_proteins[!is.na(unique_proteins)]

  if (length(unique_proteins) == 0) {
    cli::cli_alert_info("No valid protein IDs found, skipping gene symbol conversion.")
    return(var_info)
  }

  # Try to convert uniprot IDs to gene symbols
  tryCatch({
    # Use bitr to convert UNIPROT to SYMBOL, suppress warnings about mapping failures
    id_mapping <- suppressMessages(suppressWarnings(
      clusterProfiler::bitr(
        unique_proteins,
        fromType = "UNIPROT",
        toType = "SYMBOL",
        OrgDb = bitr_db
      )
    ))

    # Join with var_info to add gene symbols
    var_info <- dplyr::left_join(
      var_info,
      id_mapping,
      by = c("protein" = "UNIPROT")
    ) %>%
      dplyr::rename(gene = "SYMBOL")

    # Report conversion success
    n_converted <- sum(!is.na(var_info$gene))
    n_total <- nrow(var_info)
    perc_converted <- round(n_converted / n_total * 100, 1)
    cli::cli_alert_info(
      "Converted {.val {n_converted}} of {.val {n_total}} ({.val {perc_converted}}%) protein IDs to gene symbols."
    )
  }, error = function(e) {
    cli::cli_alert_warning(
      "Failed to convert protein IDs to gene symbols: {.val {e$message}}"
    )
  })

  var_info
}
