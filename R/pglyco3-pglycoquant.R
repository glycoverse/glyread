#' Read pGlyco3-pGlycoQuant result
#'
#' If you used pGlyco3 for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads in a pGlycoQuant result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @details
#' # Which file to use?
#'
#' You should use the "Quant.spectra.list" file in the pGlycoQuant result folder.
#' Files from pGlyco3 result folder are not needed.
#' For instructions on how to use pGlyco3 and pGlycoQuant, please refer to
#' the manual: [pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).
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
#' # Output
#'
#' This function returns a [glyexp::experiment()] object.
#'
#' The following columns could be found in the variable information tibble:
#' - `charge`: integer, charge state
#' - `peptide`: character, peptide sequence
#' - `modifications`: character, modifications other than glycosylation,
#'   separated by semicolon, e.g. `5,Carbamidomethyl[C];10,Carbamidomethyl[C]`
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures.
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `proteins`: character, protein accessions, separated by semicolon
#' - `genes`: character, gene names (symbols), separated by semicolon
#' - `protein_sites`: character, site of glycosylation on protein,
#'    separated by semicolon
#'
#' @param fp File path of the pGlycoQuant result file.
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
read_pglyco3_pglycoquant <- function(
  fp,
  sample_info = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
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
    exp <- .read_pglyco3_pglycoquant_label_free(
      fp,
      sample_info,
      glycan_type,
      sample_name_converter
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}


.read_pglyco3_pglycoquant_label_free <- function(
  fp,
  sample_info = NULL,
  glycan_type,
  sample_name_converter
) {
  # ----- Read data -----
  cli::cli_progress_step("Reading data")
  df <- .read_pglyco3_file_into_tibble(fp)
  expr_mat <- .extract_expr_mat_from_pglycoquant(df)
  samples <- .extract_sample_names_from_pglycoquant(df, sample_name_converter)
  sample_info <- .process_sample_info(sample_info, samples, glycan_type)
  var_info <- .extract_var_info_from_pglyco3(df)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable

  # ----- Parse glycan compositions and structures -----
  cli::cli_progress_step("Parsing glycan compositions and structures")
  var_info <- dplyr::mutate(
    var_info,
    glycan_composition = .convert_pglyco3_comp(.data$glycan_composition),
    glycan_structure = glyparse::parse_pglyco_struc(.data$glycan_structure)
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


# Extract variable information from pGlyco3 result
# Glycan composition and structure are not parsed here
# A "variable" column is added to the tibble
.extract_var_info_from_pglyco3 <- function(df) {
  new_names <- c(
    peptide = "Peptide",
    proteins = "Proteins",
    genes = "Genes",
    glycan_composition = "GlycanComposition",
    glycan_structure = "PlausibleStruct",
    peptide_site = "GlySite",
    protein_sites = "ProSites",
    charge = "Charge",
    modifications = "Mod"
  )
  df %>%
    dplyr::select(all_of(new_names)) %>%
    dplyr::mutate(
      modifications = stringr::str_remove(.data$modifications, ";$"),
      modifications = dplyr::if_else(is.na(.data$modifications), "", .data$modifications),
      genes = stringr::str_remove(.data$genes, ";$"),
      proteins = stringr::str_replace_all(.data$proteins, "sp\\|(\\w+)\\|\\w+", "\\1")
    ) %>%
    dplyr::mutate(variable = stringr::str_c("PSM", dplyr::row_number()), .before = 1)
}


.read_pglyco3_file_into_tibble <- function(fp) {
  col_types <- readr::cols(
    Peptide = readr::col_character(),
    Proteins = readr::col_character(),
    Genes = readr::col_character(),
    GlycanComposition = readr::col_character(),
    PlausibleStruct = readr::col_character(),
    GlySite = readr::col_integer(),
    ProSites = readr::col_character(),
    Charge = readr::col_integer(),
    Mod = readr::col_character()
  )

  # TODO: check column existence
  suppressWarnings(
    suppressMessages(readr::read_tsv(fp, col_types = col_types, progress = FALSE)),
    classes = "vroom_mismatched_column_name"
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
