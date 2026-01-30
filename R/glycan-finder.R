.parse_glycan_finder_peptide <- function(peptide) {
  stringr::str_remove_all(peptide, "\\(\\+\\d+\\.\\d+\\)")
}

.parse_glycan_finder_protein <- function(protein) {
  stringr::str_split_i(protein, stringr::fixed("|"), 1L)
}

.parse_glycan_finder_composition <- function(x) {
  # Extract monosaccharide counts from format like "(HexNAc)4(Hex)5(NeuAc)1"
  extract_n_mono <- function(comp, mono) {
    n <- stringr::str_extract(comp, paste0("\\(", mono, "\\)\\(?(\\d+)"), group = 1)
    dplyr::if_else(is.na(n), 0L, as.integer(n))
  }

  unique_x <- unique(x)
  comp_df <- purrr::map_dfc(c("HexNAc", "Hex", "NeuAc", "Fuc", "dHex", "HexA", "NeuGc", "Pent", "S", "P"), ~ {
    counts <- purrr::map_int(unique_x, extract_n_mono, mono = .x)
    tibble::tibble(!!.x := counts)
  })

  # Convert each row to a glyrepr_composition object for unique values
  unique_compositions <- purrr::pmap(comp_df, function(...) {
    counts <- c(...)
    counts <- counts[counts > 0]
    if (length(counts) == 0) {
      return(glyrepr::as_glycan_composition(integer(0)))
    }
    glyrepr::as_glycan_composition(counts)
  })

  unique_compositions <- do.call(c, unique_compositions)
  compositions <- unique_compositions[match(x, unique_x)]
  compositions
}

.read_glycan_finder_df <- function(fp) {
  col_types <- readr::cols(
    `Protein Accession` = readr::col_character(),
    Peptide = readr::col_character(),
    Glycan = readr::col_character(),
    Structure = readr::col_character(),
    `Glycan Type` = readr::col_character(),
    Start = readr::col_integer(),
    End = readr::col_integer()
  )

  suppressWarnings(
    suppressMessages(readr::read_csv(fp, col_types = col_types, progress = FALSE)),
    classes = "vroom_mismatched_column_name"
  )
}

.filter_glycan_finder_by_type <- function(df, glycan_type) {
  # Map glycan_type to GlycanFinder notation
  target_type <- if (glycan_type == "N") "N-Link" else "O-Link"

  # Filter rows that contain the target type
  df %>%
    dplyr::filter(stringr::str_detect(.data$`Glycan Type`, stringr::fixed(target_type)))
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
