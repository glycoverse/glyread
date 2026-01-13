#' Read glyco-decipher output
#'
#' Glyco-Decipher is a software for glycopeptide identification and quantification.
#' This function reads in the result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @section Which file to use?:
#' You should use the "site.csv" file in the result folder.
#' This file contains "Site" and "Glycan" columns,
#' followed by quantification result for each sample.
#'
#' @inheritSection read_pglyco3 Sample information
#'
#' @section Protein inference and uncertain sites:
#' Glyco-Decipher reports uncertain sites and proteins.
#' This function automatically performs protein inference using the
#' parsimony method to find the leader proteins.
#' This ensures each glycopeptide is uniquely mapped to a single protein, gene, and glycosite.
#' Besides, uncertain sites on proteins are assigned `NA`.
#'
#' @section Variable information:
#'
#' The following columns could be found in the variable information tibble:
#' - `protein`: character, protein accession (after protein inference)
#' - `protein_site`: integer, site of glycosylation on protein (after protein inference)
#' - `gene`: character, gene name (symbol) (after protein inference)
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#'
#' @inheritParams read_pglyco3
#' @param orgdb name of the OrgDb package to use for UniProt to gene symbol conversion.
#'  Default is "org.Hs.eg.db".
#'
#' @returns An [glyexp::experiment()] object.
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
  .validate_read_args(
    fp = fp,
    file_extensions = ".csv",
    sample_info = sample_info,
    quant_method = quant_method,
    glycan_type = glycan_type,
    sample_name_converter = sample_name_converter,
    orgdb = orgdb
  )

  # ----- Read data -----
  if (quant_method == "label-free") {
    cli::cli_progress_step("Reading data")
    df <- .read_glyco_decipher_df(fp)
    tidy_df <- .tidy_glyco_decipher(df, orgdb)
    # Add placeholder columns for standardize_variable() temporarily
    tidy_df <- dplyr::mutate(
      tidy_df,
      peptide = "N",  # Placeholder for N-glycosylation
      peptide_site = 1L
    )
    exp <- .read_template(
      tidy_df,
      sample_info,
      glycan_type,
      quant_method,
      sample_name_converter,
      composition_parser = glyrepr::as_glycan_composition,
      structure_parser = NULL,
      parse_structure = FALSE
    )
    # Remove placeholder columns so they don't appear in final var_info
    exp$var_info <- dplyr::select(exp$var_info, -"peptide", -"peptide_site")
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}

.read_glyco_decipher_df <- function(fp) {
  col_types <- readr::cols(
    Site = readr::col_character(),
    Glycan = readr::col_character()
  )
  readr::read_csv(fp, progress = FALSE, show_col_types = FALSE, col_types = col_types)
}

.tidy_glyco_decipher <- function(df, orgdb) {
  # Resolve sites
  sites <- .resolve_glyco_decipher_sites(df$Site)
  df <- df %>%
    dplyr::mutate(protein = sites$protein, protein_site = sites$protein_site) %>%
    dplyr::select(-"Site")

  # Parse glycan compositions
  df <- df %>%
    dplyr::mutate(glycan_composition = .convert_glyco_decipher_comp(.data$Glycan)) %>%
    dplyr::select(-"Glycan")

  # Add gene column
  df <- .add_gene_symbols(df, orgdb)

  # Pivot longer
  # Determine which columns to exclude from pivot_longer
  # gene column might not exist if gene conversion failed
  exclude_cols <- c("protein", "protein_site", "glycan_composition")
  if ("gene" %in% colnames(df)) {
    exclude_cols <- c(exclude_cols, "gene")
  }
  df %>%
    tidyr::pivot_longer(
      -dplyr::all_of(exclude_cols),
      names_to = "sample", values_to = "value"
    )
}

#' Resolve sites for glyco-decipher output
#'
#' @param sites The "Site" column of glyco-decipher output.
#' @returns A list of character vectors of proteins and integer vectors of sites.
#' @noRd
.resolve_glyco_decipher_sites <- function(sites) {
  sites <- stringr::str_split(sites, stringr::fixed(";"))
  proteins <- purrr::map(sites, ~ stringr::str_split_i(.x, stringr::fixed("@"), 1L))
  protein_sites <- purrr::map(sites, ~ stringr::str_split_i(.x, stringr::fixed("@"), 2L))

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
    if (length(x) == 1L) return(as.integer(x))
    return(NA_integer_)
  })
  list(protein = protein, protein_site = protein_site)
}

#' Find leader proteins for uncertain proteins in glyco-decipher output
#'
#' @param proteins A list of character vectors of proteins.
#' @returns A character vector of leader proteins.
#' @noRd
.find_glyco_decipher_leader_proteins <- function(proteins) {
  unique_proteins <- purrr::map(proteins, unique)
  unique_proteins_str <- purrr::map_chr(unique_proteins, ~ stringr::str_c(.x, collapse = ";"))
  leader_protein_idx <- .infer_proteins_internal(unique_proteins_str)
  purrr::map2_chr(unique_proteins, leader_protein_idx, ~ .x[.y])
}

.convert_glyco_decipher_comp <- function(comp) {
  comp <- stringr::str_remove(comp, "\\+.*")
  comp <- stringr::str_replace(comp, "Fuc", "dHex")
  comp
}