#' Read pGlyco3-pGlycoQuant result
#'
#' If you used pGlyco3 for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads the pGlycoQuant result file and returns an [glyexp::experiment()] object.
#'
#' @details
#' # Input
#'
#' You should use the "Quant.spectra.list" file in the pGlycoQuant result folder.
#' Files from pGlyco3 result folder are not needed.
#' For instructions on how to use pGlyco3 and pGlycoQuant, please refer to
#' the manual: [pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).
#'
#' The sample information file should be a `csv` file with the first column
#' named `sample`, and the rest of the columns are sample information.
#' The `sample` column must match the `RawName` column in the pGlyco3 result file,
#' although the order can be different.
#' You can put any useful information in the sample information file.
#' Recommended columns are:
#' - `group`: grouping or conditions, e.g. "control" or "tumor",
#'   required for most downstream analyses
#' - `batch`: batch information, required for batch effect correction
#' - `bio_replicate`:  Indicates the identity of biologically distinct samples.
#'   While the `sample` column is typically used as a unique identifier for each run,
#'   `bio_replicate` specifies whether the samples originate from the same biological source.
#'   In common experimental designs,
#'   multiple technical replicates are often derived from the same biological sample,
#'   making the bio_replicate column essential for distinguishing
#'   biological replicates from technical replicates.
#'
#' # Output
#'
#' The result is an [glyexp::experiment()] object.
#' If `parse_structure` is `TRUE` (by default),
#' a list of parsed `glycan_graph` objects will be stored in the experiment,
#' with the glycan structure strings as names.
#' It can be accessed by `exp$glycan_graphs` or `exp[["glycan_graphs"]]`.
#'
#' The following columns are in the variable information tibble:
#' - `charge`: integer, charge state
#' - `peptide`: character, peptide sequence
#' - `modifications`: character, modifications other than glycosylation,
#'   separated by semicolon, e.g. `5,Carbamidomethyl[C];10,Carbamidomethyl[C]`
#' - `glycan_composition`: character, glycan composition, e.g. "H5N4F1A1"
#' - `glycan_structure`: character, pGlyco-style structure strings, renamed from
#'   `PlausibleStruct` in the original file
#' - `n_hex`: integer, number of Hex
#' - `n_hexnac`: integer, number of HexNAc
#' - `n_neuac`: integer, number of Neu5Ac
#' - `n_neugc`: integer, number of Neu5Gc
#' - `n_fuc`: integer, number of Fuc (dHex)
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `proteins`: character, protein names, separated by semicolon
#' - `genes`: character, gene names, separated by semicolon
#' - `protein_sites`: character, site of glycosylation on protein,
#'    separated by semicolon
#'
#' If `describe_glycans` is `TRUE` (by default),
#' additional columns will be added to the variable information tibble:
#' - `glycan_type`: character, N-glycan type, either "complex", "highmannose",
#'   "hybrid", or "paucimannose"
#' - `bisecting`: logical, whether the glycan has a bisecting GlcNAc
#' - `n_antennae`: integer, number of antennae
#' - `n_core_fuc`: integer, number of core fucose
#' - `n_arm_fuc`: integer, number of arm fucose
#' - `n_gal`: integer, number of Gal
#' - `n_terminal_gal`: integer, number of terminal Gal
#'
#' Glycan compositions are reformatted into condensed format,
#' e.g. "H5N4A1F1" for "H(5)N(4)A(1)F(1)".
#'
#' @param fp File path of the pGlyco3 result file.
#' @param sample_info_fp File path of the sample information file.
#' @param name Name of the experiment. If not provided, a default name with
#'  current time will be used.
#' @param parse_structure Whether to parse glycan structures. Default is `FALSE`.
#'  Although the results from pGlyco3 provide structural information,
#'  many of these structures are speculative, and pGlyco3 doesn’t include
#'  any quality control for the structural speculation process.
#'  If you choose to enable this option,
#'  please interpret the analysis results with caution.
#'  For more rigorous structural information, we recommend using StrucGP.
#' @param describe_glycans Whether to describe glycan properties. Default is `TRUE`.
#'  If `parse_structure` is `FALSE`, this argument will be forced to `FALSE`.
#'
#' @returns An [glyexp::experiment()] object.
#'
#' @importFrom magrittr %>%
#' @importFrom tidyselect all_of
#' @importFrom rlang .data
#'
#' @export
read_pglyco3_pglycoquant <- function(
  fp,
  sample_info_fp = NULL,
  name = NULL,
  quant_method = c("label-free", "TMT"),
  parse_structure = FALSE,
  describe_glycans = TRUE
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".list")
  checkmate::assert(
    checkmate::check_null(sample_info_fp),
    checkmate::check_file_exists(sample_info_fp, access = "r", extension = ".csv")
  )
  checkmate::assert(
    checkmate::check_null(name),
    checkmate::check_character(name, len = 1, min.chars = 1)
  )
  quant_method <- rlang::arg_match(quant_method)
  checkmate::assert_flag(parse_structure)
  checkmate::assert_flag(describe_glycans)
  if (!parse_structure) {
    describe_glycans <- FALSE
  }
  if (is.null(name)) {
    name <- paste("exp", Sys.time())
  }

  # ----- Read data -----
  if (quant_method == "label-free") {
    .read_pglyco3_pglycoquant_label_free(
      fp, sample_info_fp, name, parse_structure, describe_glycans
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }
}


.read_pglyco3_pglycoquant_label_free <- function(
  fp,
  sample_info_fp = NULL,
  name = NULL,
  parse_structure = FALSE,
  describe_glycans = TRUE
) {
  # ----- Read data -----
  cli::cli_progress_step("Reading data.")
  col_types <- readr::cols(
    Peptide = readr::col_character(),
    Proteins = readr::col_character(),
    Genes = readr::col_character(),
    GlycanComposition = readr::col_character(),
    PlausibleStruct = readr::col_character(),
    GlySite = readr::col_integer(),
    ProSites = readr::col_character(),
    Charge = readr::col_integer(),
    Mod = readr::col_character()
  )
  new_names <- c(
    peptide = "Peptide",
    proteins = "Proteins",
    genes = "Genes",
    glycan_composition = "GlycanComposition",
    glycan_structure = "PlausibleStruct",
    peptide_site = "GlySite",
    protein_sites = "ProSites",
    charge = "Charge",
    modifications = "Mod"
  )

  # TODO: check column existence
  df <- suppressWarnings(
    readr::read_tsv(fp, col_types = col_types),
    classes = "vroom_mismatched_column_name"
  ) %>%
    dplyr::rename(all_of(new_names))

  samples <- unique(df$RawName)

  if (is.null(sample_info_fp)) {
    sample_info <- tibble::tibble(sample = samples)
    cli::cli_alert_info("No sample information file passed in. An empty tibble will be used.")
  } else {
    sample_info <- readr::read_csv(sample_info_fp)
    # Here we create an empty var_info tibble and an empty expr_mat matrix,
    # and check if the sample info is in correct format by creating an experiment.
    # This passes the checking responsibility to `experiment()`.
    local({
      fake_var_info <- tibble::tibble(variable = character(0))
      fake_expr_mat <- matrix(nrow = 0, ncol = length(samples), dimnames = list(NULL, samples))

      # This line will throw an error if sample_info is not in correct format.
      glyexp::experiment(name, fake_expr_mat, sample_info, fake_var_info)
    })
  }

  # ----- Clean some columns -----
  df <- df %>%
    dplyr::mutate(
      modifications = stringr::str_remove(.data$modifications, ";$"),
      modifications = dplyr::if_else(is.na(.data$modifications), "", .data$modifications),
      genes = stringr::str_remove(.data$genes, ";$")
    )

  # ----- Extract var_info -----
  var_info_cols <- names(new_names)
  var_info <- dplyr::select(df, all_of(var_info_cols))

  # ----- Parse glycan composition -----
  extract_n_mono <- function(comp, mono) {
    n <- stringr::str_extract(comp, paste0(mono, "\\((\\d+)\\)"), group = 1)
    dplyr::if_else(is.na(n), 0, as.integer(n))
  }

  var_info <- var_info %>%
    dplyr::mutate(
      n_hex = extract_n_mono(.data$glycan_composition, "H"),
      n_hexnac = extract_n_mono(.data$glycan_composition, "N"),
      n_neuac = extract_n_mono(.data$glycan_composition, "A"),
      n_neugc = extract_n_mono(.data$glycan_composition, "G"),
      n_fuc = extract_n_mono(.data$glycan_composition, "F"),
      .after = all_of("glycan_structure")
    ) %>%
    dplyr::mutate(glycan_composition = stringr::str_replace_all(
      .data$glycan_composition, "[\\(\\)]", "")
    )

  # ----- Parse glycan structure -----
  if (parse_structure) {
    cli::cli_progress_step("Parsing glycan structures.")
    glycan_structures <- unique(var_info$glycan_structure)
    glycan_graphs <- purrr::map(glycan_structures, glyparse::parse_pglyco_struc)
    names(glycan_graphs) <- glycan_structures
  } else {
    glycan_graphs <- NULL
  }

  # ----- Describe glycans -----
  if (describe_glycans) {
    cli::cli_progress_step("Extracting glycan properties (this can take long).")
    property_df <- glymotif::describe_n_glycans(glycan_graphs)
    var_info <- var_info %>%
      dplyr::left_join(property_df, by = c(glycan_structure = "glycan")) %>%
      dplyr::relocate(
        all_of(setdiff(colnames(property_df), "glycan")),
        .after = all_of("n_fuc")
      )
  }

  # ----- Pack Experiment -----
  cli::cli_progress_step("Packing experiment.")

  # Add a unique "variable" column
  var_info <- var_info %>%
    dplyr::mutate(
      variable = stringr::str_c(
        .data$peptide,
        .data$peptide_site,
        .data$glycan_composition,
        .data$glycan_structure,
        sep = "_"
      ),
      variable = make.unique(.data$variable, sep = "_"),
      .before = 1
    )

  # Extract the expression matrix
  expr_mat <- df %>%
    dplyr::select(tidyselect::starts_with("Intensity")) %>%
    as.matrix()
  colnames(expr_mat) <- stringr::str_extract(colnames(expr_mat), "Intensity\\((.*)\\)", group = 1)
  rownames(expr_mat) <- var_info$variable

  meta_data <- list(
    experiment_type = "glycoproteomics",
    glycan_type = "N-glycan",
    quantification_method = "label-free"
  )
  exp <- glyexp::experiment(name, expr_mat, sample_info, var_info, meta_data)
  exp$glycan_graphs <- glycan_graphs

  print(exp)
  invisible(exp)
}
