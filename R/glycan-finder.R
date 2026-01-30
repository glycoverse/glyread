.parse_glycan_finder_peptide <- function(peptide) {
  stringr::str_remove_all(peptide, "\\(\\+\\d+\\.\\d+\\)")
}


#' Read GlycanFinder result
#'
#' @export
read_glycan_finder <- function(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db",
  parse_structure = TRUE
) {
  stop("Not implemented")
}
