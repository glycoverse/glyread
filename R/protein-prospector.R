#' Read Protein Prospector PSM precursor intensities
#'
#' Import a tab-delimited Protein Prospector Search Compare glycopeptide report
#' into a [glyexp::GlycoproteomicSE()] object. Validated against the Protein
#' Prospector 6.5.0 report deposited in MassIVE MSV000096821.
#'
#' @section Which file to use?:
#' Use a tab-delimited report containing `DB Peptide`, `Mods`, `Glyco Mod`,
#' `Acc Num`, `Start AA`, `Fraction`, and `Intensity`. Filter identifications
#' using your search workflow before importing: this reader does not apply a
#' score or FDR threshold. Only single-glycan reports are supported. Comma-separated
#' protein accessions and their aligned start positions are resolved using the
#' package's parsimony inference; all accessions are retained in `protein_group`.
#'
#' @section Quantification limitations:
#' `Intensity` is the precursor intensity attached to the MS/MS spectrum in
#' the peak list. It is **not** an integrated XIC area or an LFQ abundance.
#' Following the other PSM readers, intensities are summed for records sharing
#' peptide, peptide site, protein, protein site, glycan composition, and sample;
#' repeated MS/MS sampling can therefore increase the summed value. Zero sums
#' become missing values. Charge states and non-glycan modifications are not
#' separate features. Each `Fraction` is a separate sample; fractions are not
#' automatically combined into biological replicates. Metadata records
#' `abundance_source = "PSM precursor intensity"`; `quant_method = "label-free"`
#' denotes unlabeled intensities, not integrated LFQ quantification.
#'
#' @section Sites and glycans:
#' A unique position in the glycan's `Mods` entry (for example,
#' `HexNAc2Hex4*@1`) determines the peptide site. Protein positions equal
#' `Start AA + peptide_site - 1`. Missing positions and alternatives such as
#' `@1|5` remain `NA`; no sequon-based localization is inferred. Unsupported
#' syntax, inconsistent glycan modifications, or out-of-range sites cause an
#' error. Supported composition tokens are HexNAc, Hex, Fuc, NeuAc, and NeuGc;
#' omitted counts mean one and Fuc is represented as dHex.
#'
#' @inheritSection read_pglyco3 Sample information
#' @inheritParams read_pglyco3
#' @param orgdb Name of the OrgDb package for gene symbol annotation.
#' @returns A [glyexp::GlycoproteomicSE()] with summed PSM precursor intensities
#'   and peptide, peptide_site, protein, protein_site, and glycan_composition
#'   row annotations (and gene when annotation is available).
#' @examples
#' # A CC0 subset of the deposited Search Compare report:
#' fp <- system.file("extdata", "protein-prospector.tsv", package = "glyread")
#' read_protein_prospector(fp, orgdb = "not.installed")
#' @export
read_protein_prospector <- function(
  fp,
  sample_info = NULL,
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db"
) {
  fp <- .read_file_input(fp, c(".txt", ".tsv"))
  .validate_read_args(
    sample_info = sample_info,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    orgdb = orgdb
  )
  df <- .read_protein_prospector_df(fp)
  tidy <- .tidy_protein_prospector(df)
  tidy <- .add_gene_symbols(tidy, orgdb)
  result <- .read_template(
    tidy,
    sample_info,
    glycan_type,
    "label-free",
    sample_name_converter,
    composition_parser = glyrepr::as_glycan_composition,
    parse_structure = FALSE
  )
  S4Vectors::metadata(result)$abundance_source <- "PSM precursor intensity"
  result
}

.read_protein_prospector_df <- function(fp) {
  # Native reports can repeat unrelated headers (e.g. Num Unique).
  df <- readr::read_tsv(
    fp,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  required <- c(
    "DB Peptide",
    "Mods",
    "Glyco Mod",
    "Acc Num",
    "Start AA",
    "Fraction",
    "Intensity"
  )
  if (any(vapply(required, \(x) sum(names(df) == x), integer(1)) != 1L)) {
    cli::cli_abort(
      "Protein Prospector requires one of each column: {.val {required}}."
    )
  }
  if (nrow(readr::problems(df)) || !nrow(df)) {
    cli::cli_abort("Protein Prospector report is empty or has malformed rows.")
  }
  df[, required]
}

.tidy_protein_prospector <- function(df) {
  required_values <- c("DB Peptide", "Glyco Mod", "Acc Num", "Fraction")
  if (
    any(vapply(
      df[required_values],
      \(x) any(is.na(x) | !nzchar(x)),
      logical(1)
    ))
  ) {
    cli::cli_abort(
      "Protein Prospector identifiers and glycans must not be missing."
    )
  }
  if (
    any(!grepl("^[A-Z]+$", df[["DB Peptide"]])) ||
      any(!grepl("^[A-Za-z0-9_-]+(,[A-Za-z0-9_-]+)*$", df[["Acc Num"]]))
  ) {
    cli::cli_abort(
      "Expected peptide sequences and protein accessions (comma-separated for aligned groups)."
    )
  }
  value <- suppressWarnings(as.numeric(df$Intensity))
  if (any(!is.na(df$Intensity) & (!is.finite(value) | value < 0))) {
    cli::cli_abort(
      "Intensity must contain non-negative finite numbers or missing values."
    )
  }
  proteins <- strsplit(df[["Acc Num"]], ",", fixed = TRUE)
  starts <- strsplit(
    ifelse(is.na(df[["Start AA"]]), "", df[["Start AA"]]),
    ",",
    fixed = TRUE
  )
  starts <- Map(
    function(pos, acc) {
      if (!length(pos)) {
        return(rep(NA_real_, length(acc)))
      }
      value <- suppressWarnings(as.numeric(pos))
      if (
        length(pos) != length(acc) ||
          any(
            !is.finite(value) |
              value < 1 |
              value != floor(value) |
              value > .Machine$integer.max - 10000
          )
      ) {
        cli::cli_abort(
          "Start AA must contain positive integer positions aligned with Acc Num."
        )
      }
      value
    },
    starts,
    proteins
  )
  groups <- gsub(",", ";", df[["Acc Num"]], fixed = TRUE)
  leader <- .infer_proteins_internal(groups)
  protein <- purrr::map2_chr(proteins, leader, \(x, i) x[i])
  start <- purrr::map2_dbl(starts, leader, \(x, i) x[i])
  site <- purrr::pmap_int(
    list(df$Mods, df[["Glyco Mod"]], df[["DB Peptide"]]),
    .protein_prospector_site
  )
  tibble::tibble(
    peptide = df[["DB Peptide"]],
    peptide_site = site,
    protein = protein,
    protein_group = groups,
    protein_site = as.integer(start + site - 1),
    glycan_composition = .protein_prospector_composition(df[["Glyco Mod"]]),
    sample = df$Fraction,
    value = value
  )
}

.protein_prospector_site <- function(mods, glycan, peptide) {
  if (is.na(mods) || !nzchar(mods)) {
    return(NA_integer_)
  }
  entries <- strsplit(mods, ";", fixed = TRUE)[[1]]
  glyco <- entries[grepl("^(HexNAc|Hex|Fuc|NeuAc|NeuGc)", entries)]
  prefix <- paste0(glycan, "*@")
  if (length(glyco) != 1L || !startsWith(glyco, prefix)) {
    cli::cli_abort(
      "Unsupported or inconsistent glycan entry in Mods: {.val {mods}}."
    )
  }
  position <- substring(glyco, nchar(prefix) + 1L)
  if (!nzchar(position)) {
    return(NA_integer_)
  }
  if (!grepl("^[0-9]+(\\|[0-9]+)*$", position)) {
    cli::cli_abort("Unsupported glycan site syntax: {.val {position}}.")
  }
  positions <- as.numeric(strsplit(position, "|", fixed = TRUE)[[1]])
  if (any(!is.finite(positions) | positions < 1 | positions > nchar(peptide))) {
    cli::cli_abort("Glycan sites must lie within the peptide sequence.")
  }
  if (length(positions) != 1L) {
    return(NA_integer_)
  }
  as.integer(positions)
}

.protein_prospector_composition <- function(x) {
  token <- "(HexNAc|Hex|Fuc|NeuAc|NeuGc)([0-9]*)"
  if (any(!grepl(paste0("^(", token, ")+$"), x))) {
    cli::cli_abort("Unsupported Protein Prospector glycan composition.")
  }
  vapply(
    x,
    function(comp) {
      parts <- stringr::str_match_all(comp, token)[[1]]
      if (anyDuplicated(parts[, 2])) {
        cli::cli_abort("Repeated glycan composition tokens.")
      }
      counts <- ifelse(parts[, 3] == "", "1", parts[, 3])
      paste0(
        ifelse(parts[, 2] == "Fuc", "dHex", parts[, 2]),
        "(",
        counts,
        ")",
        collapse = ""
      )
    },
    character(1),
    USE.NAMES = FALSE
  )
}
