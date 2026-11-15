#' Read Glyco-Decipher label-free quantification
#'
#' @section Which file to use?:
#' Pass the result directory or its `site.csv` for site-specific abundances.
#'
#' Each output row corresponds to one row of the site summary, in file
#' order. Rows remain separate even if they have the same inferred protein,
#' site and parsed composition. Blank abundances and zeros are returned as `NA`.
#'
#' @section Sample information:
#' The sample identifiers in `sample_info$sample` must match the abundance
#' column names after applying `sample_name_converter`, if supplied.
#'
#' @section Variable information:
#' `rowData` contains the inferred `protein`, the reported accessions in
#' `protein_group`, and `gene` when gene symbol mapping is available.
#' `reported_site` retains all candidate sites and
#' `protein_site` is `NA` if multiple positions remain on the inferred protein.
#' `reported_glycan` retains the original glycan label, including any mass shift, and
#' `glycan_composition` is the parsed base composition.
#'
#' @inheritParams read_pglyco3
#' @param fp Local result directory (reads `site.csv`), or a local path or
#'   HTTP(S) URL to a site summary CSV with `Site` and `Glycan` columns.
#' @param orgdb Name of the OrgDb package for mapping UniProt IDs to gene symbols.
#' @returns A [glyexp::GlycoproteomicSE()] object.
#' @examples
#' result_dir <- tempfile()
#' dir.create(result_dir)
#' writeLines(c("Site,Glycan,run1", "P02765@176,Hex(5)HexNAc(4),1250"),
#'            file.path(result_dir, "site.csv"))
#' read_glyco_decipher(result_dir)
#' @export
read_glyco_decipher <- function(
  fp,
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  orgdb = "org.Hs.eg.db"
) {
  # ----- Check arguments -----
  if (isTRUE(dir.exists(fp))) {
    fp <- file.path(fp, "site.csv")
  }
  fp <- .read_file_input(fp, ".csv")
  .validate_read_args(
    sample_info = sample_info,
    quant_method = quant_method,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    orgdb = orgdb
  )

  if (quant_method != "label-free") {
    rlang::abort("TMT quantification is not supported yet.")
  }

  cli::cli_progress_step("Reading data")
  df <- readr::read_csv(
    fp,
    col_types = readr::cols(
      Site = readr::col_character(),
      Glycan = readr::col_character(),
      .default = readr::col_double()
    ),
    name_repair = "minimal"
  )
  checkmate::assert_names(names(df), type = "unique")
  samples <- setdiff(names(df), c("Site", "Glycan"))
  if (!length(samples)) {
    cli::cli_abort("Expected at least one sample abundance column.")
  }
  readr::stop_for_problems(df)
  if (!nrow(df)) {
    cli::cli_abort("The Glyco-Decipher summary contains no features.")
  }

  abundance <- as.matrix(df[samples])
  abundance[abundance == 0] <- NA_real_
  if (!is.null(sample_name_converter)) {
    samples <- sample_name_converter(samples)
    checkmate::assert_character(
      samples,
      len = ncol(abundance),
      any.missing = FALSE,
      unique = TRUE
    )
    colnames(abundance) <- samples
  }
  sample_info <- .process_sample_info(sample_info, samples)

  var_info <- tibble::as_tibble(.resolve_glyco_decipher_sites(df$Site))
  var_info$reported_site <- df$Site
  var_info$reported_glycan <- df$Glycan
  var_info$glycan_composition <- glyrepr::as_glycan_composition(
    .convert_glyco_decipher_comp(df$Glycan)
  )
  var_info <- .add_gene_symbols(var_info, orgdb)
  var_info$variable <- paste0("GP", seq_len(nrow(df)))
  rownames(abundance) <- var_info$variable
  .new_glycoproteomic_se(
    abundance,
    sample_info,
    var_info,
    glycan_type,
    quant_method
  )
}

#' Resolve sites for glyco-decipher output
#'
#' @param sites The "Site" column of glyco-decipher output.
#' @returns A list of character vectors of proteins and integer vectors of sites.
#' @noRd
.resolve_glyco_decipher_sites <- function(sites) {
  sites <- stringr::str_split(sites, stringr::fixed(";"))
  proteins <- purrr::map(
    sites,
    ~ stringr::str_split_i(.x, stringr::fixed("@"), 1L)
  )
  protein_sites <- purrr::map(
    sites,
    ~ stringr::str_split_i(.x, stringr::fixed("@"), 2L)
  )

  # There are two situations:
  # 1. Unsure proteins: "P00738@211;P00739@153"
  # 2. Unsure sites: "P00734@416;P00734@420"
  # These two situations can happen simultaneously: "P00738@207;P00738@211;P00739@149;P00739@153"
  #
  # We first use the parsimony method to find leader proteins to deal with the first situation.
  # Then, for the second situation, the "protein_site" is assigned NA.

  # Find leader proteins
  leader_proteins <- .find_glyco_decipher_leader_proteins(proteins)
  # Keep only leader proteins
  leader_protein_mask <- purrr::map2(proteins, leader_proteins, ~ .x %in% .y)
  protein <- purrr::map2_chr(proteins, leader_protein_mask, ~ unique(.x[.y]))
  protein_sites <- purrr::map2(protein_sites, leader_protein_mask, ~ .x[.y])

  # Deal with uncertain sites
  protein_site <- purrr::map_int(protein_sites, function(x) {
    x <- unique(x)
    if (length(x) == 1L) {
      return(as.integer(x))
    }
    return(NA_integer_)
  })
  list(
    protein = protein,
    protein_site = protein_site,
    protein_group = purrr::map_chr(
      proteins,
      ~ paste(unique(.x), collapse = ";")
    )
  )
}

#' Find leader proteins for uncertain proteins in glyco-decipher output
#'
#' @param proteins A list of character vectors of proteins.
#' @returns A character vector of leader proteins.
#' @noRd
.find_glyco_decipher_leader_proteins <- function(proteins) {
  unique_proteins <- purrr::map(proteins, unique)
  unique_proteins_str <- purrr::map_chr(
    unique_proteins,
    ~ stringr::str_c(.x, collapse = ";")
  )
  leader_protein_idx <- .infer_proteins_internal(unique_proteins_str)
  purrr::map2_chr(unique_proteins, leader_protein_idx, ~ .x[.y])
}

.convert_glyco_decipher_comp <- function(comp) {
  comp <- stringr::str_remove(comp, "\\+.*")
  comp <- stringr::str_replace(comp, "Fuc", "dHex")
  comp
}
