#' Read the result from StrucGP
#'
#' StrucGP is a software for intact glycopeptide identification.
#' As StrucGP doesn't support quantification,
#' this function returns an [glyexp::experiment()] object with an expression matrix
#' filled with `NA` values.
#'
#' @details
#' # Variable information
#'
#' The following columns could be found in the variable information tibble:
#' - `peptide`: character, peptide sequence
#' - `protein`: character, protein accession
#' - `gene`: character, gene name (symbol)
#' - `protein_site`: integer, site of glycosylation on protein
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures (if `parse_structure = TRUE`).
#'
#' # Sample information
#'
#' The sample information file should be a `csv` file with the first column
#' named `sample`, and the rest of the columns being sample information.
#' The `sample` column must match the `file_name` column in the StrucGP result file,
#' although the order can be different.
#'
#' You can put any useful information in the sample information file.
#' Recommended columns are:
#' - `group`: grouping or conditions, e.g. "control" or "tumor",
#'   required for most downstream analyses
#' - `batch`: batch information, required for batch effect correction
#'
#' @param fp File path of the StrucGP result file.
#' @param sample_info File path of the sample information file (csv),
#'  or a sample information data.frame/tibble. If `NULL` (default),
#'  a simple sample information tibble will be created.
#' @param glycan_type Glycan type. Either "N" or "O". Default is "N".
#' @param parse_structure Logical. Whether to parse glycan structures.
#'  If `TRUE` (default), glycan structures are parsed and included in the
#'  `var_info` as `glycan_structure` column. If `FALSE`, structure parsing
#'  is skipped and the structure column is removed.
#'
#' @returns An [glyexp::experiment()] object with an expression matrix filled with `NA` values.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()],
#'   [glyrepr::glycan_structure()]
#' @export
read_strucgp <- function(fp, sample_info = NULL, glycan_type = "N", parse_structure = TRUE) {
  # ----- Check arguments -----
  .validate_read_args(
    fp = fp,
    file_extensions = c(".xlsx", ".xls"),
    sample_info = sample_info,
    glycan_type = glycan_type,
    parse_structure = parse_structure
  )

  # ----- Read and tidy data -----
  cli::cli_progress_step("Reading data")
  df <- readxl::read_excel(fp)
  tidy_df <- .tidy_strucgp(df, parse_structure)
  
  # ----- Extract unique samples -----
  samples <- unique(tidy_df$sample)
  
  # ----- Process sample information -----
  sample_info <- .process_sample_info(sample_info, samples, glycan_type)
  
  # ----- Create variable information -----
  cli::cli_progress_step("Creating variable information")
  var_info <- tidy_df %>%
    dplyr::select(-"sample") %>%
    dplyr::distinct() %>%
    dplyr::mutate(variable = paste0("GP", dplyr::row_number()), .before = 1)
  
  # ----- Create expression matrix with NA values -----
  expr_mat <- matrix(
    NA_real_,
    nrow = nrow(var_info),
    ncol = length(samples),
    dimnames = list(var_info$variable, samples)
  )
  
  # ----- Pack experiment -----
  cli::cli_progress_step("Creating experiment object")
  glyexp::experiment(
    expr_mat, sample_info, var_info,
    exp_type = "glycoproteomics",
    glycan_type = glycan_type,
    quant_method = "label-free"
  )
}

.tidy_strucgp <- function(df, parse_structure) {
  result <- df %>%
    janitor::clean_names() %>%
    dplyr::select(c("file_name", "peptide", "protein_id", "gene_name", "glycosite_position", "glycan_composition", "structure_coding")) %>%
    dplyr::rename(c(
      sample = "file_name",
      protein = "protein_id",
      gene = "gene_name",
      protein_site = "glycosite_position",
      glycan_structure = "structure_coding"
    )) %>%
    dplyr::mutate(
      protein_site = as.integer(.data$protein_site),
      glycan_composition = .parse_strucgp_comp(.data$glycan_composition)
    )
  
  # Parse structure only if requested
  if (parse_structure) {
    result <- dplyr::mutate(
      result,
      glycan_structure = .parse_strucgp_struc(.data$glycan_structure)
    )
  } else {
    # Remove glycan_structure column if not parsing
    result <- dplyr::select(result, -"glycan_structure")
  }
  
  result
}

.parse_strucgp_comp <- function(comp) {
  comp %>%
    stringr::str_remove("\\+.*") %>%
    stringr::str_replace("H", "Hex") %>%
    stringr::str_replace("N", "HexNAc") %>%
    stringr::str_replace("S", "NeuAc") %>%
    stringr::str_replace("G", "NeuGc") %>%
    stringr::str_replace("F", "dHex") %>%
    stringr::str_replace_all("(\\d+)", "\\(\\1\\)") %>%
    glyrepr::as_glycan_composition()
}

.parse_strucgp_struc <- function(struc) {
  struc %>%
    stringr::str_remove("\\+.*") %>%
    glyparse::parse_strucgp_struc()
}