#' Read pGlyco3-pGlycoQuant result
#'
#' If you used pGlyco3 for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads in a pGlycoQuant result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @details
#' # Input
#'
#' You should use the "Quant.spectra.list" file in the pGlycoQuant result folder.
#' Files from pGlyco3 result folder are not needed.
#' For instructions on how to use pGlyco3 and pGlycoQuant, please refer to
#' the manual: [pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).
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
#' - `bio_replicate`:  Indicates the identity of biologically distinct samples.
#'   While the `sample` column is typically used as a unique identifier for each run,
#'   `bio_replicate` specifies whether the samples originate from the same biological source.
#'   In common experimental designs,
#'   multiple technical replicates are often derived from the same biological sample,
#'   making the bio_replicate column essential for distinguishing
#'   biological replicates from technical replicates.
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
#' - `proteins`: character, protein names, separated by semicolon
#' - `genes`: character, gene names, separated by semicolon
#' - `protein_sites`: character, site of glycosylation on protein,
#'    separated by semicolon
#'
#' @param fp File path of the pGlyco3 result file.
#' @param sample_info File path of the sample information file,
#'  or a sample information tibble.
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
#'
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
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
  checkmate::assert(
    checkmate::check_null(sample_info),
    checkmate::check_file_exists(
      sample_info, access = "r", extension = ".csv"
    ),
    checkmate::check_tibble(sample_info)
  )
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  checkmate::assert(
    checkmate::check_null(sample_name_converter),
    checkmate::check_function(sample_name_converter)
  )

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
  var_info_cols <- c(
    "peptide", "proteins", "genes", "glycan_composition", "glycan_structure",
    "peptide_site", "protein_sites", "charge", "modifications"
  )

  # ----- Read data -----
  cli::cli_progress_step("Reading data")
  df <- .read_pglyco3_file_into_tibble(fp)

  if (!is.null(sample_name_converter)) {
    df$raw_name <- sample_name_converter(df$raw_name)
  }

  # Extract sample names from intensity columns to ensure consistency
  intensity_cols <- names(df)[stringr::str_starts(names(df), "Intensity")]
  samples <- stringr::str_extract(intensity_cols, "Intensity\\((.*)\\)", group = 1)
  if (!is.null(sample_name_converter)) {
    samples <- sample_name_converter(samples)
  }
  
  if (is.null(sample_info)) {
    sample_info <- tibble::tibble(sample = samples)
    cli::cli_alert_info(
      "No sample information file passed in. An empty tibble will be used."
    )
  } else {
    if (is.character(sample_info)) {
      sample_info <- suppressMessages(readr::read_csv(sample_info))
    }  # Otherwise, assume it's a tibble.
    # Here we create an empty var_info tibble and an empty expr_mat matrix,
    # and check if the sample info is in correct format
    # by creating an experiment.
    # This passes the checking responsibility to `experiment()`.
    local({
      fake_var_info <- tibble::tibble(variable = character(0))
      fake_expr_mat <- matrix(
        nrow = 0, ncol = length(samples),
        dimnames = list(NULL, samples)
      )

      # This line will throw an error if sample_info is not in correct format.
      glyexp::experiment(
        fake_expr_mat, sample_info, fake_var_info,
        exp_type = "glycoproteomics",
        glycan_type = glycan_type
      )
    })
  }

  # ----- Parse glycan compositions -----
  cli::cli_progress_step("Parsing glycan compositions")
  df <- dplyr::mutate(df, glycan_composition = .convert_glycan_composition(.data$glycan_composition))

  # ----- Parse glycan structures -----
  cli::cli_progress_step("Parsing glycan structures")
  df <- dplyr::mutate(df, glycan_structure = glyparse::parse_pglyco_struc(.data$glycan_structure))

  # ----- Pack Experiment -----
  cli::cli_progress_step("Packing experiment")
  # Add a unique "variable" column
  var_info <- df %>%
    dplyr::select(all_of(var_info_cols)) %>%
    dplyr::mutate(
      variable = stringr::str_c("PSM", dplyr::row_number()),
      .before = 1
    )

  # Extract the expression matrix
  expr_mat <- df %>%
    dplyr::select(tidyselect::starts_with("Intensity")) %>%
    as.matrix()
  expr_mat[expr_mat == 0] <- NA
  colnames(expr_mat) <- stringr::str_extract(
    colnames(expr_mat), "Intensity\\((.*)\\)", group = 1
  )
  if (!is.null(sample_name_converter)) {
    colnames(expr_mat) <- sample_name_converter(colnames(expr_mat))
  }
  rownames(expr_mat) <- var_info$variable

  exp <- glyexp::experiment(
    expr_mat, sample_info, var_info,
    exp_type = "glycoproteomics",
    glycan_type = glycan_type,
    quant_method = "label-free"
  )
  exp
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
  new_names <- c(
    raw_name = "RawName",
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

  # TODO: check column existence
  suppressWarnings(
    suppressMessages(readr::read_tsv(fp, col_types = col_types, progress = FALSE)),
    classes = "vroom_mismatched_column_name"
  ) %>%
    dplyr::rename(all_of(new_names)) %>%
    dplyr::mutate(
      modifications = stringr::str_remove(.data$modifications, ";$"),
      modifications = dplyr::if_else(
        is.na(.data$modifications), "", .data$modifications
        ),
      genes = stringr::str_remove(.data$genes, ";$"),
    )
}


.convert_glycan_composition <- function(x) {
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
  
  # Performance optimization: process unique values only
  unique_x <- unique(x)
  
  # Extract counts for each monosaccharide type for unique values
  comp_df <- purrr::map_dfc(names(pglyco_to_generic), ~ {
    counts <- purrr::map_int(unique_x, extract_n_mono, mono = .x)
    tibble::tibble(!!pglyco_to_generic[[.x]] := counts)
  })
  
  # Convert each row to a glyrepr_composition object for unique values
  unique_compositions <- purrr::pmap(comp_df, function(...) {
    counts <- c(...)
    # Keep only non-zero counts
    counts <- counts[counts > 0]
    if (length(counts) == 0) {
      # Handle empty composition - return empty composition
      return(glyrepr::as_glycan_composition(integer(0)))
    }
    glyrepr::as_glycan_composition(counts)
  })
  
  # Convert list to vector
  unique_compositions <- do.call(c, unique_compositions)
  
  # Map back to original vector using match
  compositions <- unique_compositions[match(x, unique_x)]
  
  compositions
}
