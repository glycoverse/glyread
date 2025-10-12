#' Read the result from StrucGP
#'
#' StrucGP is a software for intact glycopeptide identification.
#' As StrucGP doesn't support quantification,
#' this function only returns the processed variable information tibble,
#' instead of an [glyexp::experiment()] object like other `glyread` functions.
#'
#' @param fp File path of the StrucGP result file.
#'
#' @returns A tibble with the processed variable information.
#' @export
read_strucgp <- function(fp) {
  df <- readxl::read_excel(fp)
  df %>%
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
      glycan_composition = .parse_strucgp_comp(.data$glycan_composition),
      glycan_structure = .parse_strucgp_struc(.data$glycan_structure)
    ) %>%
    dplyr::distinct()
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