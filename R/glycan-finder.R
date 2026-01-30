.parse_glycan_finder_peptide <- function(peptide) {
  stringr::str_remove_all(peptide, "\\(\\+\\d+\\.\\d+\\)")
}

.extract_glycan_finder_peptide_site <- function(peptide, glycan_type) {
  # Remove modifications to get clean peptide
  clean_peptide <- .parse_glycan_finder_peptide(peptide)

  # Determine which residue to look for
  target_residue <- if (glycan_type == "N") "N" else "[ST]"

  # Find positions of target residues in clean peptide
  positions <- stringr::str_locate_all(clean_peptide, target_residue)[[1]]

  if (nrow(positions) == 0) {
    return(NA_integer_)
  }

  # For each position, check if the original peptide has a modification there
  for (i in seq_len(nrow(positions))) {
    pos <- positions[i, "start"]
    # Check if there's a modification at this position in original peptide
    # Pattern: X(+YYYY.YY) where X is the residue at position pos
    residue <- stringr::str_sub(clean_peptide, pos, pos)
    pattern <- paste0(residue, "\\(\\+\\d+\\.\\d+\\)")

    # Extract all modified positions from original peptide
    mods <- stringr::str_extract_all(peptide, "[A-Z]\\(\\+\\d+\\.\\d+\\)")[[1]]
    for (mod in mods) {
      mod_residue <- stringr::str_sub(mod, 1, 1)
      if (mod_residue == residue) {
        # Check if this modification is at the expected position
        # by comparing the substring before the modification
        expected_before <- stringr::str_sub(clean_peptide, 1, pos - 1)
        mod_pos <- stringr::str_locate(peptide, stringr::fixed(mod))[1, "start"]
        actual_before <- .parse_glycan_finder_peptide(stringr::str_sub(peptide, 1, mod_pos - 1))
        if (expected_before == actual_before) {
          return(as.integer(pos))
        }
      }
    }
  }

  # If no glycan modification found, return first target residue position
  as.integer(positions[1, "start"])
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
