.parse_glycan_finder_peptide <- function(peptide) {
  stringr::str_remove_all(peptide, "\\(\\+\\d+\\.\\d+\\)")
}

.extract_glycan_finder_peptide_site <- function(peptide, glycan_type) {
  # Determine which residue to look for
  target_residue <- if (glycan_type == "N") "N" else "[ST]"

  # Iterate through the peptide to find modifications
  # Pattern: X(+YYYY.YY) where X is an amino acid
  mod_pattern <- "[A-Z]\\(\\+\\d+\\.\\d+\\)"
  mods <- stringr::str_extract_all(peptide, mod_pattern)[[1]]

  if (length(mods) == 0) {
    return(NA_integer_)
  }

  # Track position in the clean peptide (without modifications)
  pos <- 0
  mod_idx <- 1
  i <- 1

  while (i <= stringr::str_length(peptide)) {
    char <- stringr::str_sub(peptide, i, i)

    # Check if this is the start of a modification
    if (mod_idx <= length(mods) && stringr::str_sub(peptide, i, i + stringr::str_length(mods[mod_idx]) - 1) == mods[mod_idx]) {
      # This is a modified residue
      residue <- stringr::str_sub(mods[mod_idx], 1, 1)
      pos <- pos + 1

      # Check if this is a target residue with a glycan modification
      # Glycan modifications are large (>500 Da), while small mods like carbamidomethylation (+57.02) are not
      mod_mass <- as.numeric(stringr::str_extract(mods[mod_idx], "\\d+\\.\\d+"))
      is_glycan <- mod_mass > 500

      if (is_glycan && stringr::str_detect(residue, target_residue)) {
        return(as.integer(pos))
      }

      # Skip past this modification
      i <- i + stringr::str_length(mods[mod_idx])
      mod_idx <- mod_idx + 1
    } else if (stringr::str_detect(char, "[A-Z]")) {
      # Regular amino acid
      pos <- pos + 1
      i <- i + 1
    } else {
      # Should not happen, but skip just in case
      i <- i + 1
    }
  }

  # If no glycan modification found, return first target residue position
  clean_peptide <- .parse_glycan_finder_peptide(peptide)
  positions <- stringr::str_locate_all(clean_peptide, target_residue)[[1]]

  if (nrow(positions) == 0) {
    return(NA_integer_)
  }

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
