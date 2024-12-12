#' Read pGlyco3-pGlycoQuant result
#'
#' If you used pGlyco3 for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads in a pGlycoQuant result file and returns a [glyexp::experiment()] object.
#' Both label-free and TMT quantification are supported.
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
#' named `sample`, and the rest of the columns being sample information.
#' The `sample` column must match the `RawName` column in the pGlyco3 result file,
#' although the order can be different.
#'
#' For results from TMT quantification, the sample information file should have
#' two additional columns besides `sample`: `raw_name` and `channel`.
#' `raw_name` refers to the raw file names,
#' matching the `RawName` column in the pGlyco3 result file.
#' `channel` refers to the TMT channel names.
#' For valid TMT channels, see the `ref_channel` argument.
#'
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
#' This function returns a [glyexp::experiment()] object.
#'
#' If `parse_structure` is `TRUE`,
#' a list of parsed `glycan_graph` objects will be stored in the experiment,
#' with the glycan structure strings as names.
#' It can be accessed by `exp$glycan_graphs` or `exp[["glycan_graphs"]]`.
#' Although the results from pGlyco3 provide structural information,
#' many of these structures are speculative, and pGlyco3 doesn’t include
#' any quality control for the structural speculation process.
#' If you choose to enable this option,
#' please interpret the analysis results with caution.
#' This is why `parse_structure` is `FALSE` by default for this function.
#'
#' The following columns could be found in the variable information tibble:
#' - `charge`: integer, charge state
#' - `peptide`: character, peptide sequence
#' - `modifications`: character, modifications other than glycosylation,
#'   separated by semicolon, e.g. `5,Carbamidomethyl[C];10,Carbamidomethyl[C]`
#' - `glycan_composition`: character, glycan composition, e.g. "H5N4F1A1"
#' - `glycan_structure`: character, pGlyco-style structure strings, renamed from
#'   `PlausibleStruct` in the original file
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `proteins`: character, protein names, separated by semicolon
#' - `genes`: character, gene names, separated by semicolon
#' - `protein_sites`: character, site of glycosylation on protein,
#'    separated by semicolon
#'
#' Some of them might not be present, depending on the aggregation method.
#' See `quant_aggr_method` argument for more information.
#'
#' # Aggregation
#'
#' pGlycoQuant quantifies glycopeptides at the PSM level.
#' Therefore, a glycopeptide could have multiple rows in the result file,
#' due to different charge states, modifications, or just duplicate identifications.
#' In this regard, we need to aggregate the quantification results
#' to the level we are interested in.
#' This function provides several aggregation methods,
#' specified by the `quant_aggr_method` argument.
#' The most commonly used method is "gf" (glycoform), which aggregates
#' all PSMs of the same combination of protein, protein site, and glycan composition.
#' The word "aggregate" here means "sum" for label-free quantification,
#' and "take the median" for TMT quantification.
#' For other aggregation methods, see the `quant_aggr_method` argument.
#'
#' @param fp File path of the pGlyco3 result file.
#' @param sample_info_fp File path of the sample information file.
#' @param name Name of the experiment. If not provided, a default name with
#'  current time will be used.
#' @param quant_method Quantification method. Either "label-free" or "TMT".
#' @param tmt_type TMT type. Required if `quant_method` is "TMT".
#'  Options are "TMT-6plex", "TMT-10plex", "TMT-11plex", "TMTpro-10plex",
#'  "TMTpro-16plex", "TMTpro-18plex".
#' @param ref_channel Reference channel for TMT normalization.
#'  A character scalar. Required if `quant_method` is "TMT".
#'  The reference channel must be one of the valid TMT channels.
#'  - TMT-6plex: 126, 127, 128, 129, 130, 131.
#'  - TMT-10plex: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131.
#'  - TMT-11plex: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N, 131C.
#'  - TMTpro-10plex: 126, 127N, 128N, 129N, 130N, 131N, 132N, 133N, 134N, 135N.
#'  - TMTpro-16plex: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N,
#'  131C, 132N, 132C, 133N, 133C, 134N.
#'  - TMTpro-18plex: 126, 127N, 127C, 128N, 128C, 129N, 129C, 130N, 130C, 131N,
#'  131C, 132N, 132C, 133N, 133C, 134N, 134C, 135N.
#' @param parse_structure Whether to parse glycan structures. Default is `FALSE`.
#'  It can only be set `TRUE` when `quant_aggr_method` is "gfs", "gps", or "det".
#'  In other cases, the structure column is dropped after aggregation.
#'  See the details section and `quant_aggr_method` argument for more information.
#' @param quant_aggr_method Aggregation method for quantification results.
#'  - "gf": Glycoform level. The default. Aggregates all PSMs of the same
#'  combination of protein, protein site, and glycan composition.
#'  - "gp": Glycopeptide level. Aggregates all PSMs of the same
#'  combination of peptide, peptide site, and glycan composition.
#'  - "gfs": Glycoform level with glycan structure.
#'  Aggregates all PSMs of the same combination of protein, protein site,
#'  and glycan structure.
#'  - "gps": Glycopeptide level with glycan structure.
#'  Aggregates all PSMs of the same combination of peptide, peptide site,
#'  and glycan structure.
#'  - "det": Detailed level. Aggregates all PSMs of the same combination of
#'  peptide, protein, gene, glycan composition, glycan structure, peptide site,
#'  protein site, charge, and modifications.
#' @param glycan_type Glycan type. Either "N" or "O". Default is "N".
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
  quant_method,
  tmt_type = NULL,
  ref_channel = NULL,
  parse_structure = FALSE,
  quant_aggr_method = c("gf", "gp", "gfs", "gps", "det"),
  glycan_type = c("N", "O")
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
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  if (!is.null(tmt_type)){
    checkmate::assert_choice(
      tmt_type,
      c("TMT-6plex", "TMT-10plex", "TMT-11plex",
        "TMTpro-10plex", "TMTpro-16plex", "TMTpro-18plex")
    )
  }
  checkmate::assert(
    checkmate::check_null(ref_channel),
    checkmate::check_string(ref_channel)
  )
  checkmate::assert_flag(parse_structure)
  quant_aggr_method <- rlang::arg_match(quant_aggr_method)
  glycan_type <- rlang::arg_match(glycan_type)
  if (is.null(name)) {
    name <- paste("exp", Sys.time())
  }
  if (!quant_aggr_method %in% c("gfs", "gps", "det")) {
    if (parse_structure) {
      parse_structure <- FALSE
      cli::cli_alert_warning(
        "parse_structure is not supported for quant_aggr_method: {.val {quant_aggr_method}}.",
        i = "The option will be ignored."
      )
    }
  }

  # ----- Read data -----
  var_info_cols <- switch (quant_aggr_method,
    gf = c("proteins", "genes", "protein_sites", "glycan_composition"),
    gp = c("peptide", "proteins", "genes", "peptide_site", "protein_sites", "glycan_composition"),
    gfs = c("proteins", "genes", "protein_sites", "glycan_composition", "glycan_structure"),
    gps = c("peptide", "proteins", "genes", "peptide_site", "protein_sites",
            "glycan_composition", "glycan_structure"),
    det = c("peptide", "proteins", "genes", "glycan_composition", "glycan_structure",
            "peptide_site", "protein_sites", "charge", "modifications")
  )
  if (quant_method == "label-free") {
    exp <- .read_pglyco3_pglycoquant_label_free(
      fp, sample_info_fp, name, glycan_type, var_info_cols)
  } else {
    exp <- .read_pglyco3_pglycoquant_tmt(
      fp, sample_info_fp, name, tmt_type, ref_channel, glycan_type, var_info_cols)
  }

  # ----- Parse glycan structure -----
  if (parse_structure) {
    glycan_structures <- unique(exp$var_info$glycan_structure)
    glycan_graphs <- purrr::map(glycan_structures, glyparse::parse_pglyco_struc)
    names(glycan_graphs) <- glycan_structures
    exp$glycan_graphs <- glycan_graphs
  }

  exp
}


.read_pglyco3_pglycoquant_label_free <- function(
  fp,
  sample_info_fp = NULL,
  name = NULL,
  glycan_type,
  var_info_cols
) {
  # ----- Read data -----
  df <- .read_pglyco3_file_into_tibble(fp)
  samples <- unique(df$raw_name)

  if (is.null(sample_info_fp)) {
    sample_info <- tibble::tibble(sample = samples)
    cli::cli_alert_info("No sample information file passed in. An empty tibble will be used.")
  } else {
    sample_info <- suppressMessages(readr::read_csv(sample_info_fp))
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

  # ----- Aggregate quantification -----
  # For label-free quantification, we sum the intensities of all
  # spectra quantified for the same variable.
  # A variable is defined as a unique combination of `names(new_names)`.
  df <- df %>%
    dplyr::summarize(
      dplyr::across(tidyselect::starts_with("Intensity"), sum),
      .by = all_of(var_info_cols)
    )
  var_info <- dplyr::select(df, all_of(var_info_cols))

  # ----- Pack Experiment -----
  # Add a unique "variable" column
  var_info <- .add_variable_column(var_info)

  # Extract the expression matrix
  expr_mat <- df %>%
    dplyr::select(tidyselect::starts_with("Intensity")) %>%
    as.matrix()
  expr_mat[expr_mat == 0] <- NA
  colnames(expr_mat) <- stringr::str_extract(colnames(expr_mat), "Intensity\\((.*)\\)", group = 1)
  rownames(expr_mat) <- var_info$variable

  meta_data <- list(
    experiment_type = "glycoproteomics",
    glycan_type = glycan_type,
    quantification_method = "label-free",
    structure_type = "pglyco"
  )
  exp <- glyexp::experiment(name, expr_mat, sample_info, var_info, meta_data)

  print(exp)
  invisible(exp)
}


.read_pglyco3_pglycoquant_tmt <- function(
  fp,
  sample_info_fp = NULL,
  name = NULL,
  tmt_type = NULL,
  ref_channel = NULL,
  glycan_type,
  var_info_cols
) {
  # ----- Read pGlyco result -----
  df <- .read_pglyco3_file_into_tibble(fp)

  # Samples are dispersed both in rows (different raw files)
  # and in columns (different TMT channels).
  raw_names <- unique(df$raw_name)
  channel_cols <- grep("^ReportIon \\(", colnames(df), value = TRUE)
  channel_mzs <- stringr::str_extract(channel_cols, "ReportIon \\((.*)\\)", group = 1)
  channels <- mz_to_tmt_channel(channel_mzs, tmt_type)
  # Rename report ion columns to TMT channel names
  df <- dplyr::rename(df, all_of(rlang::set_names(channel_cols, paste0("TMT", channels))))

  # ----- Normalize to reference channel -----
  # The reason for this step coming before sample information reading is that,
  # the reference channel is usually no longer needed after normalization,
  # so it's better to remove it from the data as soon as possible.
  if (!ref_channel %in% channels) {
    cli::cli_abort(c(
      "Reference channel not found in the data. ",
      i = "Detected channels: {.val {channels}}"
    ))
  }
  df <- df %>%
    dplyr::mutate(dplyr::across(
      all_of(paste0("TMT", channels)),
      ~ dplyr::if_else(. == 0, NA, .))
    ) %>%
    dplyr::filter(!is.na(.data[[paste0("TMT", ref_channel)]])) %>%
    dplyr::mutate(dplyr::across(
      all_of(paste0("TMT", channels)),
      ~ . / .data[[paste0("TMT", ref_channel)]])
    ) %>%
    dplyr::select(-all_of(paste0("TMT", ref_channel)))
  # remove reference channel from `channels`
  channels <- setdiff(channels, ref_channel)

  # ----- Read sample information -----
  default_sample_info <- expand.grid(
    raw_name = raw_names, channel = channels,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample = paste0("S", dplyr::row_number()))
  if (is.null(sample_info_fp)) {
    sample_info <- default_sample_info
    cli::cli_alert_info("No sample information file passed in. An empty tibble will be used.")
  } else {
    sample_info <- suppressMessages(readr::read_csv(sample_info_fp))
    if (!all(c("raw_name", "channel", "sample") %in% colnames(sample_info))) {
      cli::cli_abort("Sample information file must have columns 'raw_name', 'channel', and 'sample'.")
    }
    missing_pairs <- default_sample_info %>%
      dplyr::select(-all_of("sample")) %>%
      dplyr::left_join(sample_info, by = c("raw_name", "channel")) %>%
      dplyr::filter(is.na(.data$sample)) %>%
      dplyr::mutate(pair = paste0(.data$raw_name, ": ", .data$channel)) %>%
      dplyr::pull("pair")
    if (length(missing_pairs) > 0) {
      cli::cli_abort("Sample information is missing for these pairs: {.val {missing_pairs}}")
    }
    extra_samples <- default_sample_info %>%
      dplyr::rename(all_of(c("fake_sample" = "sample"))) %>%
      dplyr::right_join(sample_info, by = c("raw_name", "channel")) %>%
      dplyr::filter(is.na(.data$fake_sample)) %>%
      dplyr::pull("sample")
    if (length(extra_samples) > 0) {
      cli::cli_abort("Extra samples found in sample information: {.val {extra_samples}}")
    }
    if (length(unique(sample_info$sample)) != nrow(sample_info)) {
      cli::cli_abort("Sample names in the sample information file must be unique.")
    }
  }

  # ----- Aggregate quantification -----
  # For TMT quantification, we take the median value of all PSMs quantified
  # for a variable in each sample (raw name and channel pair).
  df <- df %>%
    dplyr::summarize(
      dplyr::across(tidyselect::starts_with("TMT"), ~ median(., na.rm = TRUE)),
      .by = all_of(c("raw_name", var_info_cols))
    )

  # ----- Pack experiment -----
  df <- df %>%
    tidyr::pivot_longer(
      tidyselect::starts_with("TMT"),
      names_to = "channel",
      values_to = "ratio",
      names_prefix = "TMT"
    ) %>%
    dplyr::left_join(
      sample_info[, c("raw_name", "channel", "sample")],
      by = c("raw_name", "channel")
    ) %>%
    dplyr::select(-all_of(c("raw_name", "channel"))) %>%
    tidyr::pivot_wider(names_from = "sample", values_from = "ratio")

  var_info <- df %>%
    dplyr::select(all_of(var_info_cols)) %>%
    .add_variable_column()

  expr_mat <- df %>%
    dplyr::select(-all_of(var_info_cols)) %>%
    as.matrix() %>%
    magrittr::set_rownames(var_info$variable)

  meta_data <- list(
    experiment_type = "glycoproteomics",
    glycan_type = glycan_type,
    quantification_method = "TMT",
    structure_type = "pglyco"
  )
  exp <- glyexp::experiment(name, expr_mat, sample_info, var_info, meta_data)

  print(exp)
  invisible(exp)
}


.read_pglyco3_file_into_tibble <- function(fp) {
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
    raw_name = "RawName",
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
  suppressWarnings(
    suppressMessages(readr::read_tsv(fp, col_types = col_types)),
    classes = "vroom_mismatched_column_name"
  ) %>%
    dplyr::rename(all_of(new_names)) %>%
    dplyr::mutate(
      modifications = stringr::str_remove(.data$modifications, ";$"),
      modifications = dplyr::if_else(is.na(.data$modifications), "", .data$modifications),
      genes = stringr::str_remove(.data$genes, ";$")
    ) %>%
    .clean_pglyco3_glycan_composition()
}


.clean_pglyco3_glycan_composition <- function(df) {
  extract_n_mono <- function(comp, mono) {
    n <- stringr::str_extract(comp, paste0(mono, "\\((\\d+)\\)"), group = 1)
    dplyr::if_else(is.na(n), 0, as.integer(n))
  }

  df %>%
    dplyr::mutate(
      nH = extract_n_mono(.data$glycan_composition, "H"),
      nN = extract_n_mono(.data$glycan_composition, "N"),
      nA = extract_n_mono(.data$glycan_composition, "A"),
      nG = extract_n_mono(.data$glycan_composition, "G"),
      nF = extract_n_mono(.data$glycan_composition, "F")
    ) %>%
    dplyr::mutate(glycan_composition = paste0(
      dplyr::if_else(.data$nH == 0, "", paste0("H", .data$nH)),
      dplyr::if_else(.data$nN == 0, "", paste0("N", .data$nN)),
      dplyr::if_else(.data$nF == 0, "", paste0("F", .data$nF)),
      dplyr::if_else(.data$nA == 0, "", paste0("A", .data$nA)),
      dplyr::if_else(.data$nG == 0, "", paste0("G", .data$nG))
    ), .keep = "unused")
}


.add_variable_column <- function(df) {
  # Add a unique "variable" column by concatenating all columns by "_"
  df %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character, .names = "..name..{.col}")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      variable = paste(dplyr::c_across(dplyr::starts_with("..name..")), collapse = "_"),
      .before = 1, .keep = "unused"
    ) %>%
    dplyr::ungroup()
}
