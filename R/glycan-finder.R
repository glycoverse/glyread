#' Read GlycanFinder result
#'
#' @description
#' Peaks GlycanFinder is software for intact glycopeptide identification.
#' This function reads the result file and returns a [glyexp::experiment()] object.
#'
#' @param fp File path of the GlycanFinder result CSV file.
#' @param sample_info File path of sample info or NULL for default.
#' @param quant_method Quantification method (currently only "label-free").
#' @param glycan_type Glycan type: "N", "O-GalNAc", "O-GlcNAc", etc.
#' @param sample_name_converter Function to convert sample names.
#' @param parse_structure Whether to parse structures.
#'
#' @returns glyexp::experiment object
#' @export
read_glycan_finder <- function(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  parse_structure = TRUE
) {
  rlang::abort("Not implemented yet")
}
