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
#' @inheritSection read_pglyco3 Sample information
#' @inheritSection read_pglyco3 Protein inference
#' @inheritSection read_pglyco3 Aggregation
#' @inheritSection read_pglyco3 Glycan structures
#'
#' @section Variable information:
#' The following columns could be found in the variable information tibble:
#' - `peptide`: character, peptide sequence
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accession (after protein inference)
#' - `protein_site`: integer, site of glycosylation on protein (after protein inference)
#' - `gene`: character, gene name (symbol) (after protein inference)
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures (if `parse_structure = TRUE`).
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
#' @param parse_structure Logical. Whether to parse glycan structures.
#'  If `TRUE`, glycan structures are parsed and included in the
#'  `var_info` as `glycan_structure` column. If `FALSE` (default), structure parsing
#'  is skipped and structure-related columns are removed.
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
  sample_name_converter = NULL,
  parse_structure = FALSE
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".list")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)
  checkmate::assert_logical(parse_structure, len = 1)

  # ----- Read data -----
  # Keep all PSM-level columns for variable info
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    df <- .read_pglyco3_df(fp)
    tidy_df <- .tidy_pglyco3_pglycoquant(df)
    exp <- .read_template(
      tidy_df,
      sample_info,
      glycan_type,
      quant_method,
      sample_name_converter,
      composition_parser = .convert_pglyco3_comp,
      structure_parser = glyparse::parse_pglyco_struc,
      parse_structure = parse_structure
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.tidy_pglyco3_pglycoquant <- function(df) {
  df %>%
    .convert_pglyco3_columns() %>%
    .infer_proteins_df() %>%
    .pivot_longer_pglycoquant()
}
