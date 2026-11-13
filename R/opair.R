#' Read MetaMorpheus O-Pair glycopeptide identifications
#'
#' Import native O-Pair `oglyco.psmtsv` or `glyco_psm.tsv` results as a
#' [glyexp::GlycoproteomicSE()] with a binary detection matrix, like
#' [read_strucgp()]. Fragment-ion intensities are not sample abundances.
#'
#' @details
#' Only target PSMs with `QValue <= max_q_value` are retained. Samples in the
#' input remain in the matrix even if none of their PSMs pass this filter.
#' Sample names are the exact `File Name` values, including any directory paths.
#'
#' Features are complete glycopeptides, grouped by inferred leader protein,
#' base sequence, full modified sequence and total plausible glycan composition.
#' Repeated spectra give a single detection, not a spectral count. Different
#' full sequences remain distinct. Localization probabilities are annotations,
#' not grouping keys; distinct annotations across PSMs are retained separated
#' by ` || ` in `rowData`.
#'
#' `glycan_composition` is the total composition of the complete glycopeptide,
#' never assigned independently to each of several glycosylation sites.
#' `protein_site` is populated only for an unambiguous single-protein,
#' single-glycan localization whose composition equals the total composition,
#' probability is 1, and localization level is `Level1` or `Level1b`.
#' Conflicting, ambiguous or multiple-glycan localizations have `NA` sites.
#' Original peptide/protein localization strings, localization levels and
#' plausible structures remain in `rowData`; structures are not parsed.
#' A leader protein is selected by the package's usual parsimony inference,
#' and `protein_group` retains all candidates. Gene names are unavailable.
#'
#' Supported composition symbols are H (Hex), N (HexNAc), A (NeuAc),
#' G (NeuGc), and F (dHex). Other/custom symbols fail explicitly rather
#' than being silently dropped. FragPipe's O-Pair `psm.tsv` is not this format.
#'
#' @param fp Local path or HTTP(S) URL of a native tab-separated result file.
#' @param sample_info Sample metadata CSV path/URL or data.frame with a `sample`
#'   column matching the exact `File Name` values. Defaults to `NULL`.
#' @param glycan_type Glycan type, default `"O-GalNAc"`.
#' @param max_q_value Maximum PSM QValue, between 0 and 1. Default 1
#'   retains all target PSMs; set 0.01 explicitly for a 1% PSM QValue filter.
#' @returns A [glyexp::GlycoproteomicSE()] with a binary abundance assay:
#'   1 means identified and 0 means not identified, not quantified absence.
#' @examples
#' # Read your native O-Pair output:
#' if (file.exists("oglyco.psmtsv")) {
#'   read_opair("oglyco.psmtsv")
#' }
#' @export
read_opair <- function(
  fp,
  sample_info = NULL,
  glycan_type = "O-GalNAc",
  max_q_value = 1
) {
  fp <- .read_file_input(fp, c(".psmtsv", ".tsv", ".txt"))
  .validate_read_args(sample_info = sample_info, glycan_type = glycan_type)
  checkmate::assert_number(max_q_value, lower = 0, upper = 1)
  df <- .read_opair_df(fp)
  required <- c(
    "File Name",
    "Protein Accession",
    "Base Sequence",
    "Full Sequence",
    "QValue",
    "Decoy/Contaminant/Target",
    "Plausible GlycanComposition",
    "Plausible Number Of Glycans",
    "GlycanLocalizationLevel",
    "Localized Glycans with Peptide Site Specific Probability",
    "Localized Glycans with Protein Site Specific Probability"
  )
  checkmate::assert_names(names(df), must.include = required)
  if (!nrow(df)) {
    cli::cli_abort("The O-Pair result contains no PSMs.")
  }
  if (any(is.na(df[["File Name"]]) | !nzchar(df[["File Name"]]))) {
    cli::cli_abort("O-Pair column 'File Name' contains missing values.")
  }
  if (any(!df[["Decoy/Contaminant/Target"]] %in% c("T", "D", "C"))) {
    cli::cli_abort(
      "O-Pair Decoy/Contaminant/Target must contain only T, D or C."
    )
  }
  samples <- unique(df[["File Name"]])
  sample_info <- .process_sample_info(sample_info, samples)
  q <- suppressWarnings(as.numeric(df$QValue))
  if (any(!is.finite(q) | q < 0 | q > 1)) {
    cli::cli_abort("O-Pair QValue must contain numbers between 0 and 1.")
  }
  df <- df[
    which(df[["Decoy/Contaminant/Target"]] == "T" & q <= max_q_value),
    ,
    drop = FALSE
  ]
  if (!nrow(df)) {
    cli::cli_abort("No target O-Pair PSMs pass the QValue filter.")
  }
  for (field in c(
    "Protein Accession",
    "Base Sequence",
    "Full Sequence"
  )) {
    if (any(is.na(df[[field]]) | !nzchar(df[[field]]))) {
      cli::cli_abort("O-Pair column {.val {field}} contains missing values.")
    }
  }
  tidy <- .tidy_opair(df)
  tidy$glycan_composition <- .parse_opair_comp(tidy$glycan_composition)
  keys <- c("protein", "peptide", "full_sequence", "glycan_composition")
  annotations <- c(
    "localization_level",
    "peptide_localization",
    "protein_localization",
    "plausible_structure",
    "potential_localizations",
    "site_probabilities"
  )
  vars <- tidy |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(
      protein_site = if (
        all(!is.na(.data$protein_site)) &&
          dplyr::n_distinct(.data$protein_site) == 1L
      ) {
        .data$protein_site[1L]
      } else {
        NA_integer_
      },
      dplyr::across(dplyr::all_of(annotations), .opair_annotations),
      .groups = "drop"
    )
  vars <- dplyr::left_join(
    vars,
    .protein_group_annotations(tidy, keys),
    by = keys
  )
  vars$gene <- NA_character_
  vars$variable <- paste0("GP", seq_len(nrow(vars)))
  detections <- dplyr::left_join(tidy, vars[, c(keys, "variable")], by = keys)
  mat <- matrix(
    0,
    nrow(vars),
    length(samples),
    dimnames = list(vars$variable, samples)
  )
  mat[cbind(
    match(detections$variable, vars$variable),
    match(detections$sample, samples)
  )] <- 1
  .new_glycoproteomic_se(mat, sample_info, vars, glycan_type, "label-free")
}

.opair_annotations <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x)) paste(x, collapse = " || ") else NA_character_
}

.parse_opair_comp <- function(x) {
  if (any(is.na(x) | !grepl("^([HNAGF][0-9]+)+$", x))) {
    cli::cli_abort(
      "Unsupported or missing O-Pair glycan composition; expected H/N/A/G/F followed by counts."
    )
  }
  map <- c(H = "Hex", N = "HexNAc", A = "NeuAc", G = "NeuGc", F = "dHex")
  strings <- vapply(
    stringr::str_extract_all(x, "[HNAGF][0-9]+"),
    function(tokens) {
      codes <- substr(tokens, 1, 1)
      if (anyDuplicated(codes)) {
        cli::cli_abort("Repeated symbols in O-Pair glycan composition.")
      }
      paste0(map[codes], "(", substring(tokens, 2), ")", collapse = "")
    },
    character(1)
  )
  glyrepr::as_glycan_composition(strings)
}

.opair_single_site <- function(df) {
  loc <- df[["Localized Glycans with Protein Site Specific Probability"]]
  parts <- stringr::str_match(loc, "^\\[([0-9]+),([HNAGF0-9]+),([0-9.]+)\\]$")
  valid <- !grepl("|", df[["Protein Accession"]], fixed = TRUE) &
    df[["Plausible Number Of Glycans"]] == "1" &
    df[["GlycanLocalizationLevel"]] %in% c("Level1", "Level1b") &
    parts[, 3] == df[["Plausible GlycanComposition"]] &
    suppressWarnings(as.numeric(parts[, 4])) == 1
  sites <- suppressWarnings(as.integer(parts[, 2]))
  sites[is.na(valid) | !valid | is.na(sites) | sites < 1L] <- NA_integer_
  sites
}

.read_opair_df <- function(fp) {
  # Some native exports have a trailing tab in the header only.
  lines <- readr::read_lines(fp)
  if (!length(lines)) {
    cli::cli_abort("The O-Pair result contains no PSMs.")
  }
  lines[1L] <- sub("\\t+$", "", lines[1L])
  readr::read_tsv(
    I(paste(lines, collapse = "\n")),
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE
  )
}

.tidy_opair <- function(df) {
  raw <- function(name) {
    if (!name %in% names(df)) {
      return(rep(NA_character_, nrow(df)))
    }
    df[[name]]
  }
  tidy <- tibble::tibble(
    sample = df[["File Name"]],
    protein_group = gsub("|", ";", df[["Protein Accession"]], fixed = TRUE),
    peptide = df[["Base Sequence"]],
    full_sequence = df[["Full Sequence"]],
    glycan_composition = df[["Plausible GlycanComposition"]],
    localization_level = raw("GlycanLocalizationLevel"),
    peptide_localization = raw(
      "Localized Glycans with Peptide Site Specific Probability"
    ),
    protein_localization = raw(
      "Localized Glycans with Protein Site Specific Probability"
    ),
    plausible_structure = raw("Plausible GlycanStructure"),
    potential_localizations = raw("All potential glycan localizations"),
    site_probabilities = raw("AllSiteSpecificLocalizationProbability")
  )
  tidy$protein_site <- .opair_single_site(df)
  mappings <- dplyr::distinct(tidy, .data$peptide, .data$protein_group)
  mappings$protein <- .infer_proteins(list(
    protein = mappings$protein_group
  ))$protein
  dplyr::left_join(tidy, mappings, by = c("peptide", "protein_group"))
}
