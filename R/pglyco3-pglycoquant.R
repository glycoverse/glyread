#' Read pGlyco3-pGlycoQuant result
#'
#' If you used pGlyco3 for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads in a pGlycoQuant result file and returns a [glyexp::GlycoproteomicSE()] object.
#' Currently only label-free quantification is supported.
#'
#' @details
#' # Which file to use?
#'
#' Use `Quant.site.list` (recommended) for the reported site-specific glycan
#' abundances, or `Quant.spectra.list` for peptide-level data. The reader
#' detects the format from the column headers.
#' Site files report protein, protein site and glycan composition directly;
#' protein inference, PSM aggregation and structure parsing apply only to spectra files.
#' Files from pGlyco3 result folder are not needed.
#' For instructions on how to use pGlyco3 and pGlycoQuant, please refer to
#' the manual: [pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).
#'
#' @inheritSection read_byonic_pglycoquant Sample information
#' @inheritSection read_pglyco3 Protein inference
#' @inheritSection read_pglyco3 Aggregation
#' @inheritSection read_pglyco3 Glycan structures
#'
#' @section Variable information:
#' Site files provide `protein`, `protein_site`, `glycan_composition`, and
#' `gene` when available. Spectra files can provide the following columns:
#' - `peptide`: character, peptide sequence
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accession (after protein inference)
#' - `protein_group`: character, all protein accessions before inference,
#'   separated by `;`, retained in `rowData`. Accessions are combined across
#'   records and samples when they map to the same glycopeptide.
#' - `protein_site`: integer, site of glycosylation on protein (after protein inference)
#' - `gene`: character, gene name (symbol) (after protein inference)
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures (if `parse_structure = TRUE`).
#'
#' @param orgdb Name of the OrgDb package for mapping site-file protein accessions
#'   to gene symbols. Spectra files supply their own gene annotations.
#' @param fp Local path or HTTP(S) URL of the result file.
#' @param sample_info Local path or HTTP(S) URL of the sample information CSV,
#'  or a sample information data.frame/tibble.
#' @param quant_method Quantification method. Either "label-free" or "TMT".
#' @param glycan_type Glycan type. One of "N", "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", or "O-Glc".
#'  Default is "N".
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
#' @returns An [glyexp::GlycoproteomicSE()] object.
#' @seealso [glyexp::GlycoproteomicSE()], [glyrepr::glycan_composition()],
#'   [glyrepr::glycan_structure()]
#' @export
read_pglyco3_pglycoquant <- function(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  parse_structure = FALSE,
  orgdb = "org.Hs.eg.db"
) {
  # ----- Check arguments -----
  fp <- .read_file_input(fp, ".list")
  .validate_read_args(
    sample_info = sample_info,
    quant_method = quant_method,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    parse_structure = parse_structure,
    orgdb = orgdb
  )

  # ----- Read data -----
  # Keep all PSM-level columns for variable info
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    headers <- .pglycoquant_site_headers(fp)
    if (!is.null(headers)) {
      return(.read_pglycoquant_site_table(
        fp,
        headers,
        sample_info,
        glycan_type,
        quant_method,
        sample_name_converter,
        orgdb
      ))
    }
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
