.pivot_longer_pglycoquant <- function(df) {
  df %>%
    tidyr::pivot_longer(
      tidyselect::starts_with("Intensity"),
      names_to = "sample",
      values_to = "value",
      names_pattern = "Intensity\\((.*)\\)"
    )
}

.pglycoquant_site_headers <- function(fp) {
  headers <- names(readr::read_tsv(
    fp,
    n_max = 0,
    col_types = readr::cols(),
    name_repair = "minimal"
  ))
  if (all(c("Protein", "Site", "Empty_Separator") %in% headers)) {
    headers
  }
}

# Site files already contain a wide abundance matrix in source row order.
.read_pglycoquant_site_table <- function(
  fp,
  headers,
  sample_info,
  glycan_type,
  quant_method,
  sample_name_converter,
  orgdb
) {
  glycan_col <- grep("^Glycan\\(.*\\)$", headers, value = TRUE)
  if (length(glycan_col) != 1L) {
    cli::cli_abort(
      "Expected one Glycan(...) column in the pGlycoQuant site file."
    )
  }
  sample_cols <- which(
    seq_along(headers) > match("Empty_Separator", headers) & nzchar(headers)
  )
  samples <- headers[sample_cols]
  if (!length(samples)) {
    cli::cli_abort("Expected sample columns after Empty_Separator.")
  }
  types <- rep("_", length(headers))
  types[headers %in% c("Protein", glycan_col)] <- "c"
  types[headers == "Site"] <- "i"
  types[sample_cols] <- "d"
  df <- readr::read_tsv(
    fp,
    col_types = paste0(types, collapse = ""),
    name_repair = "minimal"
  )
  df <- .drop_reverse_rows(df, "Protein", "pGlycoQuant")
  abundance <- as.matrix(df[, samples, drop = FALSE])
  if (!is.null(sample_name_converter)) {
    samples <- sample_name_converter(samples)
  }
  checkmate::assert_character(
    samples,
    len = length(sample_cols),
    any.missing = FALSE,
    unique = TRUE
  )
  sample_info <- .process_sample_info(sample_info, samples)
  abundance[abundance == 0] <- NA_real_
  colnames(abundance) <- samples

  var_info <- tibble::tibble(
    variable = paste0("GP", seq_len(nrow(df))),
    protein = .extract_uniprot_accession(df$Protein),
    protein_site = df$Site,
    glycan_composition = .parse_pglycoquant_site_composition(
      df[[glycan_col]],
      glycan_col
    )
  ) |>
    .add_gene_symbols(orgdb)
  rownames(abundance) <- var_info$variable
  .new_glycoproteomic_se(
    abundance,
    sample_info,
    var_info,
    glycan_type,
    quant_method
  )
}

.parse_pglycoquant_site_composition <- function(x, header) {
  codes <- trimws(strsplit(
    sub("^Glycan\\((.*)\\)$", "\\1", header),
    ",",
    fixed = TRUE
  )[[1]])
  if (
    anyDuplicated(codes) ||
      !all(codes %in% c("H", "N", "A", "G", "F", "aH", "pH"))
  ) {
    cli::cli_abort("Unsupported monosaccharide labels in {.val {header}}.")
  }
  # Capture one count per header label, in the order used by this export.
  pattern <- paste0(
    "^\\[\\s*",
    paste(rep("(\\d+)", length(codes)), collapse = "\\s+"),
    "\\s*\\]$"
  )
  distinct <- unique(x)
  counts <- stringr::str_match(distinct, pattern)
  if (anyNA(counts)) {
    cli::cli_abort(
      "Expected one nonnegative count per monosaccharide in {.val {header}}."
    )
  }
  pglyco <- apply(counts[, -1, drop = FALSE], 1, function(row) {
    paste0(codes, "(", row, ")", collapse = "")
  })
  .convert_pglyco3_comp(pglyco)[match(x, distinct)]
}
