#' Read GlycanFinder result
#'
#' @description
#' GlycanFinder is a software for intact glycopeptide identification.
#' This function reads in the result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @section Which file to use?:
#' You should use the `lfq/lfq.protein-glycopeptides.csv` file from the GlycanFinder
#' result folder. This file contains the quantification information for glycopeptides.
#'
#' @section Sample information:
#' The sample information file should be a `csv` file with the first column
#' named `sample`, and the rest of the columns being sample information.
#' The `sample` column must match the sample names in the Area columns of the
#' GlycanFinder result file (e.g., "C_3", "H_3"), although the order can be different.
#'
#' @section Variable information:
#' The following columns could be found in the variable information tibble:
#' - `peptide`: character, peptide sequence
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accession
#' - `protein_site`: integer, site of glycosylation on protein
#' - `gene`: character, gene name (symbol)
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures (if `parse_structure = TRUE`).
#'
#' @param fp File path of the GlycanFinder result file.
#' @param sample_info File path of the sample information file (csv),
#'  or a sample information data.frame/tibble.
#' @param quant_method Quantification method. Either "label-free" or "TMT".
#' @param glycan_type Glycan type. One of "N", "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", or "O-Glc".
#'  Default is "N".
#' @param sample_name_converter A function to convert sample names from file paths.
#'  The function should take a character vector of old sample names
#'  and return new sample names.
#'  Note that sample names in `sample_info` should match the new names.
#'  If NULL, original names are kept.
#' @param orgdb Name of the OrgDb package to use for UniProt to gene symbol conversion.
#'  Default is "org.Hs.eg.db".
#' @param parse_structure Logical. Whether to parse glycan structures.
#'  If `TRUE` (default), glycan structures are parsed and included in the
#'  `var_info` as `glycan_structure` column. If `FALSE`, structure parsing
#'  is skipped and structure-related columns are removed.
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()],
#'   [glyrepr::glycan_structure()]
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
  # ----- Check arguments -----
  .validate_read_args(
    fp = fp,
    file_extensions = ".csv",
    sample_info = sample_info,
    quant_method = quant_method,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    parse_structure = parse_structure,
    orgdb = orgdb
  )

  # ----- Read data -----
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    df <- .read_glycan_finder_df(fp)
    tidy_df <- .tidy_glycan_finder(df, glycan_type, orgdb)
    exp <- .read_template(
      tidy_df,
      sample_info,
      glycan_type,
      quant_method,
      sample_name_converter,
      composition_parser = .parse_glycan_finder_composition,
      structure_parser = .parse_glycan_finder_structure,
      parse_structure = parse_structure
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.parse_glycan_finder_peptide <- function(peptide) {
  stringr::str_remove_all(peptide, "\\(\\+\\d+\\.\\d+\\)")
}

.parse_glycan_finder_protein <- function(protein) {
  stringr::str_split_i(protein, stringr::fixed("|"), 1L)
}

.parse_glycan_finder_structure <- function(x) {
  # Wrapper around glyparse::parse_pglyco_struc that handles semicolon-separated structures
  # Some GlycanFinder structures contain multiple structures separated by ';'
  # We parse only the first structure in such cases
  parse_single <- function(struc) {
    if (is.na(struc) || struc == "") {
      return(NA)
    }
    # If structure contains semicolon, take only the first part
    if (stringr::str_detect(struc, ";")) {
      struc <- stringr::str_split_i(struc, ";", 1L)
    }
    # Parse using glyparse
    tryCatch(
      glyparse::parse_pglyco_struc(struc),
      error = function(e) NA
    )
  }

  unique_x <- unique(x)
  unique_structures <- purrr::map(unique_x, parse_single)
  unique_structures <- do.call(c, unique_structures)
  structures <- unique_structures[match(x, unique_x)]
  structures
}

.parse_glycan_finder_composition <- function(x) {
  # Extract monosaccharide counts from format like "(HexNAc)4(Hex)5(NeuAc)1"
  # Note: We map Fuc -> dHex and Pent -> Pen to match glyrepr's naming conventions
  extract_n_mono <- function(comp, mono) {
    n <- stringr::str_extract(comp, paste0("\\(", mono, "\\)\\(?(\\d+)"), group = 1)
    dplyr::if_else(is.na(n), 0L, as.integer(n))
  }

  unique_x <- unique(x)
  # Use dHex instead of Fuc to avoid mixing generic and concrete monosaccharide types
  # Use Pen instead of Pent (glyrepr uses "Pen" for pentose)
  # S and P are substituents, not monosaccharides - we skip them for now
  mono_names <- c("HexNAc", "Hex", "NeuAc", "dHex", "HexA", "NeuGc", "Pen")
  search_names <- c("HexNAc", "Hex", "NeuAc", "Fuc", "HexA", "NeuGc", "Pent")

  comp_df <- purrr::map2_dfc(search_names, mono_names, ~ {
    counts <- purrr::map_int(unique_x, extract_n_mono, mono = .x)
    tibble::tibble(!!.y := counts)
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

.tidy_glycan_finder <- function(df, glycan_type, orgdb = "org.Hs.eg.db") {
  df %>%
    .filter_glycan_finder_by_type(glycan_type) %>%
    .convert_glycan_finder_columns(glycan_type) %>%
    .pivot_glycan_finder_long() %>%
    .add_gene_symbols(orgdb)
}

.convert_glycan_finder_columns <- function(df, glycan_type) {
  df %>%
    dplyr::mutate(
      peptide = .parse_glycan_finder_peptide(.data$Peptide),
      protein = .parse_glycan_finder_protein(.data$`Protein Accession`),
      peptide_site = purrr::map_int(.data$Peptide, .extract_glycan_finder_peptide_site, glycan_type = glycan_type),
      protein_site = .data$peptide_site + .data$Start - 1L,
      glycan_composition = .select_glycan_element(.data$`Glycan Type`, .data$Glycan, glycan_type),
      glycan_structure = .select_glycan_element(.data$`Glycan Type`, .data$Structure, glycan_type)
    )
}

# Helper function to select the correct glycan element from semicolon-separated values
# For mixed glycan types (e.g., "N-Link;O-Link"), selects the element matching the target type
.select_glycan_element <- function(glycan_types, values, target_glycan_type) {
  # Determine which index to select (1 for N, 2 for O)
  target_index <- if (target_glycan_type == "N") 1L else 2L

  purrr::map2_chr(glycan_types, values, function(gt, val) {
    # Check if this is a mixed type row
    if (stringr::str_detect(gt, ";")) {
      # Split by semicolon and select the appropriate element
      elements <- stringr::str_split_1(val, ";")
      if (length(elements) >= target_index) {
        return(elements[target_index])
      }
    }
    # Return as-is for non-mixed types or if split fails
    val
  })
}

.pivot_glycan_finder_long <- function(df) {
  # Find Area columns (sample abundance columns)
  area_cols <- colnames(df)[stringr::str_detect(colnames(df), "^Area ")]
  sample_names <- stringr::str_remove(area_cols, "^Area ")

  # Pivot to long format
  df %>%
    dplyr::select(
      "peptide", "peptide_site", "protein", "protein_site",
      "glycan_composition", "glycan_structure",
      dplyr::all_of(area_cols)
    ) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(area_cols),
      names_to = "sample",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      sample = stringr::str_remove(.data$sample, "^Area "),
      value = dplyr::if_else(.data$value == 0, NA_real_, .data$value)
    )
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
