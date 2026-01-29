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

.tidy_glycan_finder <- function(df, glycan_type, parse_structure) {
  # Clean names
  tidy_df <- janitor::clean_names(df)

  # Extract sample columns (Area *)
  area_cols <- grep("^area_", names(tidy_df), value = TRUE)

  # Pivot to long format
  long_df <- tidy_df %>%
    dplyr::select(dplyr::all_of(c(
      "glycan", "structure", "glycan_type", "peptide",
      "protein_accession", area_cols
    ))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(area_cols),
      names_to = "sample_raw",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      sample = stringr::str_remove(.data$sample_raw, "^area_"),
      # Replace 0 with NA
      value = dplyr::if_else(.data$value == 0, NA_real_, .data$value)
    ) %>%
    dplyr::select(-"sample_raw")

  # Filter by glycan type
  long_df <- .filter_glycan_type(long_df, glycan_type)

  # Rename columns
  long_df <- long_df %>%
    dplyr::rename(
      glycan_composition = "glycan",
      glycan_structure = "structure"
    )

  # Parse structures if requested
  if (parse_structure) {
    cli::cli_progress_step("Parsing glycan structures")
    long_df <- dplyr::mutate(
      long_df,
      glycan_structure = glyparse::parse_strucgp_struc(.data$glycan_structure)
    )
  } else {
    long_df <- dplyr::select(long_df, -"glycan_structure")
  }

  long_df
}

.filter_glycan_type <- function(df, glycan_type) {
  # Map glycan_type parameter to GlycanFinder format
  glycan_type_map <- c(
    "N" = "N-Link",
    "O-GalNAc" = "O-Link",
    "O-GlcNAc" = "O-Link",
    "O-Man" = "O-Link",
    "O-Fuc" = "O-Link",
    "O-Glc" = "O-Link"
  )

  target_type <- glycan_type_map[[glycan_type]]
  if (is.null(target_type)) {
    cli::cli_abort("Unknown glycan type: {.val {glycan_type}}")
  }

  # Filter rows where glycan_type column contains target type
  dplyr::filter(df, stringr::str_detect(.data$glycan_type, target_type))
}
