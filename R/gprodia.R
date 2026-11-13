#' Read GproDIA results
#'
#' @description
#' Read a GproDIA long-format TSV report into a
#' [glyexp::GlycoproteomicSE()] object. Only label-free quantification is supported.
#'
#' @details
#' Use a `legacy_merged` or `legacy_split` report from GproDIA, or its
#' TRIC-aligned `*_aligned.tsv` output, not a matrix export. Required columns are
#' `filename`, `run_id`, `transition_group_id`, `Sequence`, `GlycanSite`,
#' `GlycanComposition`, `ProteinName`, `ProteinGlycoSite`, `Intensity`,
#' `m_score`, `decoy`, `decoy_peptide`, and `decoy_glycan`.
#'
#' Quantification uses `Intensity` (MS2 fragment peak area), as in GproDIA's
#' matrix export, not `aggr_prec_Peak_Area` (MS1 precursor area). Decoys are
#' removed. When multiple peak groups occur for one run and transition group,
#' the row with the lowest `m_score` is retained, matching GproDIA's matrix
#' export. Ties use the first row. No additional score threshold is applied:
#' perform confidence filtering in GproDIA/TRIC before importing.
#'
#' Protein accessions and sites separated by `/` or `;` are paired by position
#' and resolved using glyread's parsimony protein inference. Original accessions
#' are retained in `protein_group`. Gene names are not inferred. The peptide
#' sequence uses `N` in place of GproDIA's glycosylated-asparagine marker `J`.
#' Charge states, peptide modifications and glycan structural isomers are
#' summed at peptide, peptide-site, protein-site and glycan-composition level.
#' Glycan structures are not imported. Zero totals become missing values.
#'
#' Sample names are taken from `filename`, including its extension.
#' `sample_info` must contain a `sample` column matching those names (or the
#' converted names when `sample_name_converter` is supplied).
#'
#' @inheritParams read_pglyco3
#' @returns A [glyexp::GlycoproteomicSE()] object.
#' @export
read_gprodia <- function(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL
) {
  fp <- .read_file_input(fp, ".tsv")
  .validate_read_args(
    sample_info,
    quant_method,
    glycan_type,
    sample_name_converter
  )
  if (quant_method != "label-free") {
    rlang::abort("TMT quantification is not supported for GproDIA.")
  }
  df <- .read_gprodia_df(fp)
  tidy_df <- .tidy_gprodia(df)
  .read_template(
    tidy_df,
    sample_info,
    glycan_type,
    quant_method,
    sample_name_converter,
    composition_parser = .convert_pglyco3_comp,
    parse_structure = FALSE
  )
}

.read_gprodia_df <- function(fp) {
  readr::read_tsv(
    fp,
    col_types = readr::cols(.default = readr::col_character()),
    progress = FALSE
  )
}

.tidy_gprodia <- function(df) {
  required <- c(
    "filename",
    "run_id",
    "transition_group_id",
    "Sequence",
    "GlycanSite",
    "GlycanComposition",
    "ProteinName",
    "ProteinGlycoSite",
    "Intensity",
    "m_score",
    "decoy",
    "decoy_peptide",
    "decoy_glycan"
  )
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    cli::cli_abort("GproDIA report is missing columns: {.val {missing}}.")
  }
  for (column in c("decoy", "decoy_peptide", "decoy_glycan")) {
    if (anyNA(df[[column]]) || any(!df[[column]] %in% c("0", "1"))) {
      cli::cli_abort("GproDIA {.val {column}} must contain only 0 or 1.")
    }
  }
  df <- dplyr::filter(
    df,
    .data$decoy == "0",
    .data$decoy_peptide == "0",
    .data$decoy_glycan == "0"
  )
  if (!nrow(df)) {
    cli::cli_abort("GproDIA report contains no target rows.")
  }
  if (
    anyNA(df[required]) ||
      any(vapply(df[required], \(x) any(!nzchar(x)), logical(1)))
  ) {
    cli::cli_abort("Required GproDIA fields must not be missing or empty.")
  }
  run_files <- dplyr::distinct(df, .data$run_id, .data$filename)
  if (anyDuplicated(run_files$run_id)) {
    cli::cli_abort(
      "Each GproDIA run_id must map to exactly one filename; import one export report at a time."
    )
  }
  for (column in c("Intensity", "m_score", "GlycanSite")) {
    value <- suppressWarnings(as.numeric(df[[column]]))
    if (any(!is.finite(value)) || (column != "m_score" && any(value < 0))) {
      cli::cli_abort(
        "GproDIA {.val {column}} contains invalid numeric values."
      )
    }
    df[[column]] <- value
  }
  if (
    any(
      df$GlycanSite < 1 |
        df$GlycanSite != floor(df$GlycanSite) |
        df$GlycanSite > nchar(df$Sequence)
    )
  ) {
    cli::cli_abort(
      "GproDIA GlycanSite must be a valid one-based peptide position."
    )
  }
  if (any(!grepl("^((H|N|A|G|F)\\([0-9]+\\))+$", df$GlycanComposition))) {
    cli::cli_abort("Unsupported GproDIA glycan composition notation.")
  }
  df <- df |>
    dplyr::group_by(.data$run_id, .data$transition_group_id) |>
    dplyr::slice_min(.data$m_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  proteins <- strsplit(
    gsub("/", ";", df$ProteinName, fixed = TRUE),
    ";",
    fixed = TRUE
  )
  sites <- strsplit(
    gsub("/", ";", df$ProteinGlycoSite, fixed = TRUE),
    ";",
    fixed = TRUE
  )
  if (
    any(lengths(proteins) != lengths(sites)) ||
      any(!grepl("^[1-9][0-9]*$", unlist(sites))) ||
      any(!nzchar(unlist(proteins)))
  ) {
    cli::cli_abort(
      "GproDIA protein and site lists must be non-empty and paired by position."
    )
  }
  proteins <- lapply(proteins, \(x) sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", x))
  tidy <- tibble::tibble(
    peptide = gsub("J", "N", df$Sequence, fixed = TRUE),
    peptide_site = as.integer(df$GlycanSite),
    proteins = vapply(proteins, paste, character(1), collapse = ";"),
    protein_sites = vapply(sites, paste, character(1), collapse = ";"),
    genes = vapply(
      proteins,
      \(x) paste(rep("NA", length(x)), collapse = ";"),
      character(1)
    ),
    glycan_composition = df$GlycanComposition,
    sample = df$filename,
    value = df$Intensity
  )
  tidy <- .infer_proteins_df(tidy)
  tidy$gene <- NULL
  tidy
}
