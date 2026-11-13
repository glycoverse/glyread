#' Read GlycReSoft glycopeptide identification CSV exports
#'
#' Read label-free N-glycopeptide exports in the format verified against the
#' GlycReSoft 0.4.23 paper's public search results (doi:10.6084/m9.figshare.24578857).
#' This does not read spectrum-match CSVs or SQLite databases. Other software
#' versions have not been validated; incompatible development exports with
#' `index:sequon` site annotations are rejected. Each `protein_name` must be a
#' UniProt-style FASTA title (`sp|accession|...` or `tr|accession|...`).
#'
#' A positive `group_id`, scoped to its file and `analysis`, identifies one
#' shared chromatographic signal. Its `total_signal` is counted once, not once
#' per alternative identification. No additional q-value filtering is applied.
#' `total_signal` measures MS1 chromatographic signal, not PSM counts.
#' Zero signal and features not observed in a sample are missing abundance
#' (`NA`). Apply desired q-value/FDR filtering upstream.
#' Different groups with the same feature
#' identity are summed. Feature identity includes the inferred leader protein,
#' all candidate modified glycopeptides/compositions, and the leader's candidate
#' peptide positions. It excludes sample-local group IDs and localization
#' probabilities, allowing equivalent candidate sets to align across samples.
#'
#' Protein inference uses parsimony. Occupied sites come from the
#' `N(N-Glycosylation)` markers, not all potential sequons in
#' `n_glycosylation_sites`. Protein coordinates are converted from zero-based
#' starts to one-based sites. Ambiguous or multiple occupied sites are `NA`.
#' A site is reported only when all retained candidates agree and their
#' localization probabilities are one (within numerical tolerance `1e-8`).
#' Missing or lower probabilities retain the signal with an unknown site.
#' Multiple candidate glycan compositions produce one feature with an `NA`
#' composition, retaining its abundance and all alternatives.
#'
#' @param fp Character vector of local identification CSV paths or HTTP(S) URLs.
#'   Sample names come from `analysis`, which must not occur in multiple input
#'   files. A file can contain several analyses.
#' @inheritParams read_msfragger
#' @param sample_name_converter Optional function transforming the `analysis`
#'   sample names into unique names matching `sample_info`.
#' @returns A [glyexp::GlycoproteomicSE()] with label-free abundance. Row metadata
#'   includes `feature_signature`, `protein_group`, `candidate_glycopeptides`,
#'   `candidate_compositions`, and `source_annotations` (a list of original
#'   source rows with file and sample identifiers). Candidate annotations remain
#'   available even when composition or occupied site is unknown.
#' @examples
#' # Small excerpt of the authors' original search results:
#' fp <- system.file("extdata", "glycresoft.csv", package = "glyread")
#' read_glycresoft(fp)
#' @export
read_glycresoft <- function(
  fp,
  sample_info = NULL,
  sample_name_converter = NULL
) {
  checkmate::assert_character(
    fp,
    min.len = 1L,
    any.missing = FALSE,
    unique = TRUE
  )
  .validate_read_args(
    sample_info,
    sample_name_converter = sample_name_converter
  )
  files <- vapply(fp, .read_file_input, character(1), extensions = ".csv")
  input <- lapply(seq_along(files), function(i) {
    x <- readr::read_csv(
      files[[i]],
      col_types = readr::cols(.default = "c"),
      name_repair = "check_unique",
      show_col_types = FALSE
    )
    if (nrow(readr::problems(x))) {
      cli::cli_abort("Malformed CSV: parser reported inconsistent records.")
    }
    .validate_glycresoft(x)
    x$source_file <- fp[[i]]
    x$source_sample <- x$analysis
    x
  })
  data <- dplyr::bind_rows(input)
  mapping <- unique(data[c("analysis", "source_file")])
  if (anyDuplicated(mapping$analysis)) {
    cli::cli_abort(
      "An analysis occurs in multiple files; supply each sample only once."
    )
  }
  groups <- dplyr::group_rows(dplyr::group_by(
    data,
    .data$source_file,
    .data$analysis,
    .data$group_id
  ))
  entities <- lapply(groups, function(i) .glycresoft_entity(data[i, ]))
  # Infer once per distinct candidate mapping, not once per replicate.
  mappings <- vapply(
    entities,
    function(x) paste(x$proteins, collapse = ";"),
    ""
  )
  evidence <- vapply(
    entities,
    function(x) paste(sort(unique(x$peptides)), collapse = ";"),
    ""
  )
  unique_map <- !duplicated(paste(evidence, mappings, sep = "\r"))
  choices <- .infer_proteins_internal(mappings[unique_map])
  leaders <- mapply(
    function(x, i) strsplit(x, ";", fixed = TRUE)[[1]][i],
    mappings[unique_map],
    choices,
    USE.NAMES = FALSE
  )
  leaders <- leaders[match(
    paste(evidence, mappings, sep = "\r"),
    paste(evidence[unique_map], mappings[unique_map], sep = "\r")
  )]
  entities <- Map(.glycresoft_feature, entities, leaders)
  samples <- unique(data$analysis)
  converted <- if (is.null(sample_name_converter)) {
    samples
  } else {
    sample_name_converter(samples)
  }
  checkmate::assert_character(
    converted,
    len = length(samples),
    any.missing = FALSE,
    unique = TRUE,
    min.chars = 1L
  )
  sample_info <- .process_sample_info(sample_info, converted)
  signatures <- vapply(entities, `[[`, "", "feature_signature")
  indices <- split(seq_along(entities), signatures)
  keys <- sort(names(indices))
  abundance <- matrix(
    NA_real_,
    length(keys),
    length(samples),
    dimnames = list(paste0("GP", seq_along(keys)), converted)
  )
  features <- lapply(seq_along(keys), function(k) {
    items <- entities[indices[[keys[[k]]]]]
    for (s in samples) {
      values <- vapply(
        items[vapply(items, function(x) identical(x$sample, s), TRUE)],
        `[[`,
        0,
        "value"
      )
      if (length(values)) {
        abundance[k, match(s, samples)] <<- if (sum(values) == 0) {
          NA_real_
        } else {
          sum(values)
        }
      }
    }
    first <- items[[1]]
    rows <- dplyr::bind_rows(lapply(items, `[[`, "rows"))
    tibble::tibble(
      variable = paste0("GP", k),
      peptide = first$peptide,
      peptide_site = if (
        length(unique(vapply(items, `[[`, 0L, "peptide_site"))) == 1L
      ) {
        first$peptide_site
      } else {
        NA_integer_
      },
      protein = first$protein,
      protein_site = if (
        length(unique(vapply(items, `[[`, 0L, "protein_site"))) == 1L
      ) {
        first$protein_site
      } else {
        NA_integer_
      },
      gene = first$gene,
      glycan_composition = first$composition,
      feature_signature = keys[[k]],
      protein_group = paste(
        sort(unique(unlist(lapply(items, `[[`, "proteins")))),
        collapse = ";"
      ),
      candidate_glycopeptides = list(sort(unique(rows$glycopeptide))),
      candidate_compositions = list(first$compositions),
      source_annotations = list(rows)
    )
  })
  var_info <- dplyr::bind_rows(features)
  var_info$glycan_composition <- glyrepr::as_glycan_composition(
    var_info$glycan_composition
  )
  result <- .new_glycoproteomic_se(
    abundance,
    sample_info,
    var_info,
    "N",
    "label-free"
  )
  S4Vectors::metadata(result)$abundance_source <- "GlycReSoft total_signal"
  result
}

.validate_glycresoft <- function(x) {
  required <- c(
    "glycopeptide",
    "analysis",
    "total_signal",
    "start_time",
    "end_time",
    "apex_time",
    "peptide_start",
    "peptide_end",
    "protein_name",
    "n_glycosylation_sites",
    "mass_shifts",
    "group_id",
    "localizations"
  )
  if (!all(required %in% names(x)) || !nrow(x)) {
    cli::cli_abort(
      "Expected a nonempty GlycReSoft identification CSV with group_id and localization annotations."
    )
  }
  if (
    anyNA(x[c("glycopeptide", "analysis", "protein_name", "group_id")]) ||
      any(!nzchar(x$analysis))
  ) {
    cli::cli_abort("Missing GlycReSoft identity fields.")
  }
  for (column in c("group_id", "peptide_start", "peptide_end")) {
    if (anyNA(x[[column]]) || any(!grepl("^[0-9]+$", x[[column]]))) {
      cli::cli_abort(
        "{column} must contain nonnegative integer identifiers/coordinates."
      )
    }
  }
  if (any(as.numeric(x$group_id) < 1)) {
    cli::cli_abort(
      "group_id must be positive; missing or fallback IDs cannot identify shared signals."
    )
  }
  value <- suppressWarnings(as.numeric(x$total_signal))
  if (anyNA(value) || any(!is.finite(value) | value < 0)) {
    cli::cli_abort("total_signal must be finite and nonnegative.")
  }
  site <- x$n_glycosylation_sites
  if (any(!is.na(site) & !grepl("^[0-9]+(;[0-9]+)*$", site))) {
    cli::cli_abort(
      "Unsupported n_glycosylation_sites format; expected zero-based integer lists, not index:sequon annotations."
    )
  }
}

.glycresoft_parse <- function(x) {
  parts <- stringr::str_match(x, "^(.+)\\{([^{}]+)\\}$")
  if (anyNA(parts)) {
    cli::cli_abort("Malformed GlycReSoft glycopeptide annotation.")
  }
  annotated <- parts[, 2]
  peptide <- gsub("\\([^()]*\\)", "", annotated)
  if (any(!grepl("^[A-Z]+$", peptide))) {
    cli::cli_abort("Unsupported modified peptide syntax.")
  }
  sites <- lapply(annotated, function(s) {
    tokens <- stringr::str_extract_all(s, "[A-Z](?:\\([^()]*\\))*")[[1]]
    if (!identical(paste0(tokens, collapse = ""), s)) {
      cli::cli_abort("Unsupported modified peptide syntax.")
    }
    hit <- which(grepl("(N-Glycosylation)", tokens, fixed = TRUE))
    if (!length(hit) || any(substr(tokens[hit], 1, 1) != "N")) {
      cli::cli_abort("Expected occupied N-glycosylation markers.")
    }
    hit
  })
  compositions <- vapply(
    parts[, 3],
    function(s) {
      tokens <- strsplit(s, ";", fixed = TRUE)[[1]]
      fields <- stringr::str_match(
        trimws(tokens),
        "^([A-Za-z][A-Za-z0-9]*):([0-9]+)$"
      )
      if (anyNA(fields) || anyDuplicated(fields[, 2])) {
        cli::cli_abort("Malformed GlycReSoft glycan composition.")
      }
      aliases <- c(
        Fuc = "dHex",
        Hex = "Hex",
        HexNAc = "HexNAc",
        Neu5Ac = "NeuAc",
        Neu5Gc = "NeuGc"
      )
      names <- unname(aliases[fields[, 2]])
      if (anyNA(names)) {
        cli::cli_abort(
          "Unsupported GlycReSoft monosaccharide; refusing to discard a component."
        )
      }
      as.character(glyrepr::as_glycan_composition(paste0(
        names,
        "(",
        fields[, 3],
        ")",
        collapse = ""
      )))
    },
    ""
  )
  list(
    peptide = peptide,
    sites = sites,
    compositions = unname(compositions),
    annotated = annotated
  )
}

.glycresoft_entity <- function(rows) {
  for (column in c(
    "total_signal",
    "start_time",
    "end_time",
    "apex_time",
    "mass_shifts"
  )) {
    if (length(unique(rows[[column]])) != 1L) {
      cli::cli_abort(
        "Conflicting {column} within a GlycReSoft shared-signal group."
      )
    }
  }
  parsed <- .glycresoft_parse(rows$glycopeptide)
  starts <- as.numeric(rows$peptide_start)
  ends <- as.numeric(rows$peptide_end)
  if (
    any(ends - starts != nchar(parsed$peptide)) ||
      any(ends > .Machine$integer.max)
  ) {
    cli::cli_abort("Peptide coordinates do not match the annotated sequence.")
  }
  for (i in seq_len(nrow(rows))) {
    raw_sites <- rows$n_glycosylation_sites[[i]]
    if (!is.na(raw_sites)) {
      potential <- as.numeric(strsplit(raw_sites, ";", fixed = TRUE)[[1]])
      if (any(potential < starts[[i]] | potential >= ends[[i]])) {
        cli::cli_abort(
          "Potential glycosylation sites fall outside peptide coordinates."
        )
      }
    }
  }
  proteins <- stringr::str_match(
    rows$protein_name,
    "^>?(?:sp|tr)\\|([^|;[:space:]]+)\\|"
  )[, 2]
  if (anyNA(proteins)) {
    cli::cli_abort("Expected UniProt protein_name annotations.")
  }
  genes <- stringr::str_match(rows$protein_name, "(?:^| )GN=([^ ]+)")[, 2]
  list(
    rows = rows,
    sample = rows$analysis[[1]],
    value = as.numeric(rows$total_signal[[1]]),
    proteins = sort(unique(proteins)),
    row_proteins = proteins,
    genes = genes,
    peptides = parsed$peptide,
    sites = parsed$sites,
    starts = starts,
    annotated = parsed$annotated,
    row_compositions = parsed$compositions,
    compositions = sort(unique(parsed$compositions))
  )
}

.glycresoft_feature <- function(entity, leader) {
  e <- entity
  selected <- which(e$row_proteins == leader)
  one <- function(x) if (length(unique(x)) == 1L) unique(x) else NA
  tuples <- paste(e$annotated, e$row_compositions, sep = "{")
  positions <- paste(e$annotated[selected], e$starts[selected], sep = "@")
  e$feature_signature <- paste(
    leader,
    paste(sort(unique(tuples)), collapse = "|"),
    paste(sort(unique(positions)), collapse = "|"),
    sep = "::"
  )
  e$protein <- leader
  e$peptide <- as.character(one(e$peptides))
  e$gene <- as.character(one(e$genes[selected]))
  peptide_sites <- e$sites[selected]
  e$peptide_site <- if (all(lengths(peptide_sites) == 1L)) {
    as.integer(one(unlist(peptide_sites)))
  } else {
    NA_integer_
  }
  protein_sites <- Map(
    function(s, start) s + start,
    peptide_sites,
    e$starts[selected]
  )
  e$protein_site <- if (all(lengths(protein_sites) == 1L)) {
    as.integer(one(unlist(protein_sites)))
  } else {
    NA_integer_
  }
  localized <- vapply(
    selected,
    function(i) {
      annotation <- e$rows$localizations[[i]]
      if (is.na(annotation) || !nzchar(annotation)) {
        return(FALSE)
      }
      annotations <- strsplit(annotation, ";", fixed = TRUE)[[1]]
      annotations <- annotations[startsWith(annotations, "N-Glycosylation:")]
      if (!length(annotations)) {
        return(FALSE)
      }
      fields <- stringr::str_match(
        annotations,
        "^N-Glycosylation:([0-9]+):([0-9.eE+-]+)$"
      )
      if (anyNA(fields)) {
        cli::cli_abort("Unsupported localization annotation.")
      }
      probability <- suppressWarnings(as.numeric(fields[, 3]))
      if (
        anyNA(probability) ||
          any(
            !is.finite(probability) | probability < 0 | probability > 1 + 1e-8
          )
      ) {
        cli::cli_abort("Invalid localization probability.")
      }
      occupied <- e$sites[[i]] - 1L
      match_index <- match(occupied, as.integer(fields[, 2]))
      !anyNA(match_index) && all(abs(probability[match_index] - 1) <= 1e-8)
    },
    TRUE
  )
  if (!all(localized)) {
    e$peptide_site <- NA_integer_
    e$protein_site <- NA_integer_
  }
  e$composition <- if (length(e$compositions) == 1L) {
    e$compositions
  } else {
    NA_character_
  }
  e
}
