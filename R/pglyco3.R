#' Read pGlyco3 Result
#'
#' This function reads pGlyco3 quantification result,
#' and returns an [glyexp::experiment()] object.
#'
#' @details
#' # Input
#'
#' You should use the `txt` file in pGlyco3 result folder with "Quant" in the
#' file name. It is something like `pGlycoDB-GP-FDR-Pro-Quant.txt`.
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
#' - `glycan_composition`: character, glycan composition, e.g. "H(5)N(4)A(1)F(1)"
#' - `glycan_structure`: character, pGlyco-style structure strings, renamed from
#'   `PlausibleStruct` in the original file
#' - `n_hex`: integer, number of Hex
#' - `n_hexnac`: integer, number of HexNAc
#' - `n_neuac`: integer, number of NeuAc
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
#' Besides, glycan compositions are reformatted into condensed format,
#' e.g. "H5N4F1S1" for "H(5)N(4)A(1)F(1)".
#' Note that "A" and "G" in glycan structures are replaced by "S"
#' if `differ_a_g` is `FALSE`.
#'
#' @param fp File path of the pGlyco3 result file.
#' @param sample_info_fp File path of the sample information file.
#' @param name Name of the experiment. If not provided, a default name with
#'  current time will be used.
#' @param quantify_on Quantify on "mono" or "sum". If "mono", the "MonoArea"
#'  column will be used for quantification. If "sum", the "IsotopeArea" column
#'  will be used for quantification. Default is "mono".
#' @param differ_a_g Whether to differentiate "A" (Neu5Ac) and "G" (Neu5Gc).
#'  If `FALSE`, "A" and "G" will be replaced by "S" in compositions and structures.
#'  Default is `FALSE`.
#' @param parse_structure Whether to parse glycan structures. Default is `TRUE`.
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
read_pglyco3 <- function(
  fp,
  sample_info_fp = NULL,
  name = NULL,
  quantify_on = c("mono", "sum"),
  differ_a_g = FALSE,
  parse_structure = TRUE,
  describe_glycans = TRUE
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".txt")
  checkmate::assert(
    checkmate::check_null(sample_info_fp),
    checkmate::check_file_exists(sample_info_fp, access = "r", extension = ".csv")
  )
  checkmate::assert(
    checkmate::check_null(name),
    checkmate::check_character(name, len = 1, min.chars = 1)
  )
  quantify_on <- rlang::arg_match(quantify_on)
  checkmate::assert_flag(parse_structure)
  checkmate::assert_flag(describe_glycans)
  if (!parse_structure) {
    describe_glycans <- FALSE
  }
  if (is.null(name)) {
    name <- paste("exp", Sys.time())
  }

  # ----- Read data -----
  cli::cli_progress_step("Reading data.")
  col_types <- readr::cols_only(
    RawName = readr::col_character(),
    Peptide = readr::col_character(),
    Proteins = readr::col_character(),
    Genes = readr::col_character(),
    `Glycan(H,N,A,F)` = readr::col_character(),
    PlausibleStruct = readr::col_character(),
    GlySite = readr::col_integer(),
    ProSites = readr::col_character(),
    Charge = readr::col_integer(),
    Mod = readr::col_character(),
    MonoArea = readr::col_double(),
    IsotopeArea = readr::col_double()
  )
  new_names <- c(
    sample = "raw_name",
    glycan_structure = "plausible_struct",
    peptide_site = "gly_site",
    protein_sites = "pro_sites",
    modifications = "mod"
  )
  df <- readr::read_delim(fp, delim = "\t", col_types = col_types) %>%
    janitor::clean_names() %>%
    dplyr::rename(all_of(new_names))

  if (is.null(sample_info_fp)) {
    sample_info <- tibble::tibble(sample = unique(df$sample))
    cli::cli_alert_info("No sample information file passed in. An empty tibble will be used.")
  } else {
    sample_info <- readr::read_csv(sample_info_fp)
  }

  # ----- Check sample info -----
  # Here we create an empty var_info tibble and an empty expr_mat matrix,
  # and check if the sample info is in correct format by creating an experiment.
  # This passes the checking responsibility to `experiment()`.
  local({
    fake_var_info <- tibble::tibble(variable = character(0))
    samples <- unique(df$sample)
    fake_expr_mat <- matrix(nrow = 0, ncol = length(samples), dimnames = list(NULL, samples))

    # This line will throw an error if sample_info is not in correct format.
    glyexp::experiment(name, fake_expr_mat, sample_info, fake_var_info)
  })

  # ----- Deal with A and G -----
  if (!differ_a_g) {
    df <- dplyr::mutate(
      df, glycan_structure = stringr::str_replace_all(.data$glycan_structure, "A", "S")
    )
    # `sia_sym` will be used in the "Parse glycan composition" section
    sia_sym <- "S"
  } else {
    sia_sym <- "A"
  }

  # ----- Decide Quantification Column -----
  if (quantify_on == "mono") {
    df <- df %>%
      dplyr::rename(all_of(c(area = "mono_area"))) %>%
      dplyr::select(-all_of("isotope_area"))
  } else {
    df <- df %>%
      dplyr::rename(all_of(c(area = "isotope_area"))) %>%
      dplyr::select(-all_of("mono_area"))
  }

  # ----- Clean some columns -----
  df <- df %>%
    dplyr::mutate(
      modifications = stringr::str_remove(.data$modifications, ";$"),
      modifications = dplyr::if_else(is.na(.data$modifications), "", .data$modifications),
      genes = stringr::str_remove(.data$genes, ";$")
    )

  # ----- Parse glycan composition -----
  df <- df %>%
    tidyr::separate_wider_delim(
      cols = all_of("glycan_h_n_a_f"),
      delim = " ",
      names = c("n_hex", "n_hexnac", "n_neuac", "n_fuc")
    ) %>%
    dplyr::mutate(dplyr::across(
      all_of(c("n_hex", "n_hexnac", "n_neuac", "n_fuc")), as.integer
    )) %>%
    dplyr::relocate(
      all_of(c("n_hex", "n_hexnac", "n_neuac", "n_fuc")),
      .after = all_of("glycan_structure")
    ) %>%
    dplyr::mutate(glycan_composition = stringr::str_c(
      dplyr::if_else(.data$n_hex > 0, stringr::str_c("H", .data$n_hex), ""),
      dplyr::if_else(.data$n_hexnac > 0, stringr::str_c("N", .data$n_hexnac), ""),
      dplyr::if_else(.data$n_fuc > 0, stringr::str_c("F", .data$n_fuc), ""),
      dplyr::if_else(.data$n_neuac > 0, stringr::str_c(sia_sym, .data$n_neuac), "")
    ), .before = dplyr::all_of("glycan_structure"))

  # ----- Parse glycan structure -----
  if (parse_structure) {
    cli::cli_progress_step("Parsing glycan structures.")
    glycan_structures <- unique(df$glycan_structure)
    glycan_graphs <- purrr::map(glycan_structures, glyparse::parse_pglyco_struc)
    names(glycan_graphs) <- glycan_structures
  } else {
    glycan_graphs <- NULL
  }

  # ----- Describe glycans -----
  if (describe_glycans) {
    cli::cli_progress_step("Extracting glycan properties (this can take long).")
    property_df <- glymotif::describe_n_glycans(glycan_graphs)
    df <- df %>%
      dplyr::left_join(property_df, by = c(glycan_structure = "glycan")) %>%
      dplyr::relocate(
        all_of(setdiff(colnames(property_df), "glycan")),
        .after = all_of("n_fuc")
      )
  }

  # ----- Pack Experiment -----
  cli::cli_progress_step("Packing experiment.")
  df_wide <- df %>%
    tidyr::pivot_wider(
      names_from = all_of("sample"),
      values_from = all_of("area"),
      values_fn = sum
    ) %>%
    dplyr::mutate(
      variable = stringr::str_c(
        .data$peptide,
        .data$peptide_site,
        .data$glycan_structure,
        sep = "_"
      ),
      variable = make.unique(.data$variable, sep = "_")
    )
  var_info <- df_wide %>%
    dplyr::select(all_of(c("variable", setdiff(colnames(df), c("sample", "area")))))
  expr_mat <- df_wide %>%
    dplyr::select(all_of(c("variable", unique(df$sample)))) %>%
    tibble::column_to_rownames("variable") %>%
    as.matrix()
  exp <- glyexp::experiment(name, expr_mat, sample_info, var_info)
  exp$glycan_graphs <- glycan_graphs

  print(exp)
  invisible(exp)
}
