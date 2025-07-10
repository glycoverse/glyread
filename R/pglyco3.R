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
    exp <- .read_pglyco3_label_free(
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

.read_pglyco3_label_free <- function(
  fp,
  sample_info = NULL,
  glycan_type,
  sample_name_converter
) {
  # ----- Read data -----
  cli::cli_progress_step("Reading data")
  df <- .read_pglyco3_file_into_tibble(fp) %>%
    .convert_pglyco3_columns()
  # Apply sample name converter to the data frame if provided
  if (!is.null(sample_name_converter)) {
    df <- dplyr::mutate(df, RawName = sample_name_converter(.data$RawName))
  }
  sample_names <- unique(df$RawName)

  sample_info <- .process_sample_info(sample_info, sample_names, glycan_type)

  # ----- Aggregate PSMs to glycopeptides -----
  cli::cli_progress_step("Aggregating PSMs to glycopeptides")
  aggregated_result <- df %>%
    dplyr::rename(all_of(c(sample = "RawName", value = "MonoArea"))) %>%
    .aggregate_long()
  var_info <- aggregated_result$var_info
  expr_mat <- aggregated_result$expr_mat

  # ----- Protein inference -----
  cli::cli_progress_step("Performing protein inference")
  var_info <- .infer_proteins(var_info)

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
