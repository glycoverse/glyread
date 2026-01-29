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
  # Validate arguments
  .validate_read_args(
    fp = fp,
    file_extensions = ".csv",
    sample_info = sample_info,
    quant_method = quant_method,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    parse_structure = parse_structure
  )

  # Read data
  cli::cli_progress_step("Reading GlycanFinder data")
  df <- suppressMessages(readr::read_csv(fp, progress = FALSE))

  # Tidy data
  tidy_df <- .tidy_glycan_finder(df, glycan_type, parse_structure)

  # Apply sample name converter
  if (!is.null(sample_name_converter)) {
    tidy_df <- dplyr::mutate(tidy_df, sample = sample_name_converter(.data$sample))
  }

  # Get unique samples
  samples <- unique(tidy_df$sample)

  # Process sample info
  sample_info_out <- .process_sample_info(sample_info, samples, glycan_type)

  # Create variable info (before parsing compositions to preserve join keys)
  cli::cli_progress_step("Creating variable information")
  var_info <- tidy_df %>%
    dplyr::select(-"sample", -"value") %>%
    dplyr::distinct() %>%
    dplyr::mutate(variable = paste0("GP", dplyr::row_number()), .before = 1)

  # Join tidy_df with var_info to get variable IDs for expression matrix
  # Exclude columns that are in tidy_df but not in var_info (like protein_site and end)
  join_by <- setdiff(names(tidy_df), c("sample", "value", "protein_site", "end"))
  # Conditionally select glycan_structure only if it exists
  select_cols <- c("variable", "glycan_composition", "peptide", "protein", "glycan_type")
  if ("glycan_structure" %in% names(var_info)) {
    select_cols <- c(select_cols, "glycan_structure")
  }
  tidy_df_with_var <- dplyr::left_join(
    tidy_df,
    dplyr::select(var_info, dplyr::all_of(select_cols)),
    by = join_by
  )

  # Parse glycan compositions
  var_info <- dplyr::mutate(
    var_info,
    glycan_composition = .parse_glycan_finder_comp(.data$glycan_composition)
  )

  # Pivot to wide format for expression matrix
  expr_mat <- matrix(
    NA_real_,
    nrow = nrow(var_info),
    ncol = length(samples),
    dimnames = list(var_info$variable, samples)
  )

  # Fill expression matrix
  for (i in seq_len(nrow(tidy_df_with_var))) {
    var_id <- tidy_df_with_var$variable[i]
    sample_id <- tidy_df_with_var$sample[i]
    val <- tidy_df_with_var$value[i]
    expr_mat[var_id, sample_id] <- val
  }

  # Pack experiment
  cli::cli_progress_step("Creating experiment object")
  exp <- glyexp::experiment(
    expr_mat, sample_info_out, var_info,
    exp_type = "glycoproteomics",
    glycan_type = glycan_type,
    quant_method = quant_method
  )

  # Standardize variable IDs
  glyexp::standardize_variable(exp)
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
      "protein_accession", "start", "end", area_cols
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
      glycan_structure = "structure",
      protein = "protein_accession",
      protein_site = "start"
    )

  # Clean protein accession (remove isoform info after "|")
  long_df <- dplyr::mutate(
    long_df,
    protein = stringr::str_split(.data$protein, "[|]") %>%
      purrr::map_chr(1)
  )

  # Parse structures if requested
  if (parse_structure) {
    cli::cli_progress_step("Parsing glycan structures")
    # Take first alternative if multiple structures exist (separated by ";")
    long_df <- dplyr::mutate(
      long_df,
      glycan_structure = purrr::map_chr(
        stringr::str_split(.data$glycan_structure, ";"),
        1
      ),
      glycan_structure = glyparse::auto_parse(.data$glycan_structure)
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

.parse_glycan_finder_comp <- function(comp) {
  # Take first alternative if multiple compositions exist (separated by ";")
  comp <- purrr::map_chr(stringr::str_split(comp, ";"), 1)

  # Convert from "(HexNAc)4(Hex)5(NeuAc)1" to "HexNAc(4)Hex(5)NeuAc(1)"
  # Also convert "Fuc" to "dHex" for consistency (glyrepr expects consistent naming)
  comp %>%
    stringr::str_remove_all("[()]") %>%
    stringr::str_replace_all("([A-Za-z]+)([0-9]+)", "\\1(\\2)") %>%
    stringr::str_replace("Fuc", "dHex") %>%
    glyrepr::as_glycan_composition()
}
