#' Read pGlyco3 result
#'
#' pGlyco3 is a software for intact glycopeptide identification and quantification.
#' This function reads in the result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @details
#' # Which file to use?
#'
#' You should use the result file from pGlyco3 that contains quantification information.
#' The file should have columns including `RawName`, `MonoArea`, `Peptide`, `Proteins`,
#' `Genes`, `GlycanComposition`, `PlausibleStruct`, `GlySite`, and `ProSites`.
#'
#' # Sample information
#'
#' The sample information file should be a `csv` file with the first column
#' named `sample`, and the rest of the columns being sample information.
#' The `sample` column must match the `RawName` column in the pGlyco3 result file,
#' although the order can be different.
#'
#' You can put any useful information in the sample information file.
#' Recommended columns are:
#' - `group`: grouping or conditions, e.g. "control" or "tumor",
#'   required for most downstream analyses
#' - `batch`: batch information, required for batch effect correction
#'
#' # Protein inference
#'
#' By default, this function automatically performs protein inference using the
#' parsimony method to resolve shared glycopeptides. This converts the plural
#' columns (`proteins`, `genes`, `protein_sites`) to singular equivalents
#' (`protein`, `gene`, `protein_site`).
#'
#' # Aggregation
#'
#' pGlyco3 performs quantification on the PSM level.
#' This level of information is too detailed for most downstream analyses.
#' This function aggregate PSMs into glycopeptides through summation.
#' For each glycopeptide (unique combination of "peptide", "peptide_site", "protein", "protein_site",
#' "gene", "glycan_composition", "glycan_structure"),
#' we sum up the quantifications of all PSMs that belong to this glycopeptide.
#'
#' # Output
#'
#' This function returns a [glyexp::experiment()] object.
#'
#' The following columns could be found in the variable information tibble:
#' - `peptide`: character, peptide sequence
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accession (after protein inference)
#' - `protein_site`: integer, site of glycosylation on protein (after protein inference)
#' - `gene`: character, gene name (symbol) (after protein inference)
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures.
#'
#' @param fp File path of the pGlyco3 result file.
#' @param sample_info File path of the sample information file (csv),
#'  or a sample information data.frame/tibble.
#' @param quant_method Quantification method. Either "label-free" or "TMT".
#' @param glycan_type Glycan type. Either "N" or "O". Default is "N".
#' @param sample_name_converter A function to convert sample names from file paths.
#'  The function should take a character vector of old sample names
#'  and return new sample names.
#'  Note that sample names in `sample_info` should match the new names.
#'  If NULL, original names are kept.
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()],
#'   [glyrepr::glycan_structure()]
#' @export
read_pglyco3 <- function(
  fp,
  sample_info = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".txt")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)

  # ----- Read data -----
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    df <- .read_pglyco3_df(fp)
    tidy_df <- .tidy_pglyco3(df)
    exp <- .read_template(
      tidy_df,
      sample_info,
      glycan_type,
      quant_method,
      sample_name_converter,
      composition_parser = .convert_pglyco3_comp,
      structure_parser = glyparse::parse_pglyco_struc
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.tidy_pglyco3 <- function(df) {
  df %>%
    .convert_pglyco3_columns() %>%
    .infer_proteins_df() %>%
    dplyr::rename(sample = "RawName", value = "MonoArea")
}

.read_pglyco3_df <- function(fp) {
  col_types <- readr::cols(
    Peptide = readr::col_character(),
    Proteins = readr::col_character(),
    Genes = readr::col_character(),
    GlycanComposition = readr::col_character(),
    PlausibleStruct = readr::col_character(),
    GlySite = readr::col_integer(),
    ProSites = readr::col_character(),
  )

  suppressWarnings(
    suppressMessages(readr::read_tsv(fp, col_types = col_types, progress = FALSE)),
    classes = "vroom_mismatched_column_name"
  )
}

.convert_pglyco3_columns <- function(df) {
  new_names <- c(
    peptide = "Peptide",
    proteins = "Proteins",
    genes = "Genes",
    glycan_composition = "GlycanComposition",
    glycan_structure = "PlausibleStruct",
    peptide_site = "GlySite",
    protein_sites = "ProSites"
  )
  df %>%
    dplyr::rename(all_of(new_names)) %>%
    dplyr::mutate(
      genes = stringr::str_remove(.data$genes, ";$"),
      proteins = stringr::str_replace_all(.data$proteins, "sp\\|(\\w+)\\|\\w+", "\\1")
    )
}

.convert_pglyco3_comp <- function(x) {
  # Define mapping from pGlyco3 notation to generic monosaccharides
  pglyco_to_generic <- c(
    "H" = "Hex",      # Hexose
    "N" = "HexNAc",   # N-Acetylhexosamine  
    "A" = "NeuAc",    # N-Acetylneuraminic acid
    "G" = "HexA",     # Hexuronic acid
    "F" = "dHex"      # Deoxyhexose (Fucose)
  )
  
  extract_n_mono <- function(comp, mono) {
    n <- stringr::str_extract(comp, paste0(mono, "\\((\\d+)\\)"), group = 1)
    dplyr::if_else(is.na(n), 0L, as.integer(n))
  }

  unique_x <- unique(x)
  comp_df <- purrr::map_dfc(names(pglyco_to_generic), ~ {
    counts <- purrr::map_int(unique_x, extract_n_mono, mono = .x)
    tibble::tibble(!!pglyco_to_generic[[.x]] := counts)
  })
  
  # Convert each row to a glyrepr_composition object for unique values
  unique_compositions <- purrr::pmap(comp_df, function(...) {
    counts <- c(...)
    counts <- counts[counts > 0]
    if (length(counts) == 0) {
      return(glyrepr::as_glycan_composition(integer(0)))
    }
    glyrepr::as_glycan_composition(counts)
  })

  unique_compositions <- do.call(c, unique_compositions)
  compositions <- unique_compositions[match(x, unique_x)]
  compositions
}
