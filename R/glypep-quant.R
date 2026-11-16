#' Read GlyPep-Quant label-free quantification results
#'
#' Import a GlyPep-Quant CSV produced from Glyco-Decipher or pGlyco3
#' identifications. The column headers identify the summary format:
#'
#' - `site.csv` reports a site-specific glycan per row in `Site` and `Glycan`.
#' - `glycopeptide.csv` reports a glycopeptide per row in `Peptide`, `Glycan`,
#'   and `Protein`.
#'
#' All remaining columns contain sample abundances. Empty cells and zeros
#' become `NA`. Rows are retained in file order, including distinct modified
#' forms of the same peptide or glycan.
#'
#' @seealso [Published GlyPep-Quant result files (PXD057799)](https://www.ebi.ac.uk/pride/archive/projects/PXD057799)
#' @inheritParams read_glyco_decipher
#' @param quant_method Quantification method. Only `"label-free"` is supported.
#' @returns A [glyexp::GlycoproteomicSE()] with one row per input record.
#'   Row metadata include the inferred leader `protein`, the candidate
#'   accessions in `protein_group`, and `protein_site` (an integer, or `NA`
#'   for absent or ambiguous coordinates). `reported_glycan` retains the
#'   original label, including any mass shift, alongside `glycan_composition`.
#'   Site tables retain `reported_site`. Glycopeptide tables retain the
#'   original `reported_protein` and `reported_peptide` labels and provide the
#'   unmodified sequence in `peptide`.
#' @examples
#' fp <- tempfile(fileext = ".csv")
#' writeLines(c("Peptide,Glycan,Protein,HeLa_01",
#'              "AAFHNAVAVVIYNNK,Hex(8)HexNAc(2),Q86XS8@135;,831800"), fp)
#' se <- read_glypep_quant(fp)
#' unlink(fp)
#' @export
read_glypep_quant <- function(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db"
) {
  .read_glyco_summary(
    fp,
    sample_info,
    quant_method,
    glycan_type,
    sample_name_converter,
    orgdb,
    col_types = NULL
  )
}

.glypep_quant_col_types <- function(fp) {
  headers <- names(readr::read_csv(
    fp,
    n_max = 0,
    col_types = readr::cols(),
    name_repair = "minimal"
  ))
  if (all(c("Peptide", "Glycan", "Protein") %in% headers)) {
    return(readr::cols(
      Peptide = "c",
      Glycan = "c",
      Protein = "c",
      .default = "d"
    ))
  }
  if (all(c("Site", "Glycan") %in% headers)) {
    return(readr::cols(Site = "c", Glycan = "c", .default = "d"))
  }
  cli::cli_abort(
    "Expected GlyPep-Quant columns Site and Glycan, or Peptide, Glycan and Protein."
  )
}
