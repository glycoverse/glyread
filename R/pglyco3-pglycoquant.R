#' Read pGlyco3-pGlycoQuant result
#'
#' If you used pGlyco3 for intact glycopeptide identification,
#' and used pGlycoQuant for quantification, this is the function for you.
#' It reads in a pGlycoQuant result file and returns a [glyexp::experiment()] object.
#' Currently only label-free quantification is supported.
#'
#' @details
#' # Which file to use?
#'
#' You should use the "Quant.spectra.list" file in the pGlycoQuant result folder.
#' Files from pGlyco3 result folder are not needed.
#' For instructions on how to use pGlyco3 and pGlycoQuant, please refer to
#' the manual: [pGlycoQuant](https://github.com/Power-Quant/pGlycoQuant/blob/main/Manual%20for%20pGlycoQuant_v202211.pdf).
#'
#' # Sample information
#'
#' The sample information file should be a `csv` file with the first column
#' named `sample`, and the rest of the columns being sample information.
#' The `sample` column must match the `RawName` column in the pGlyco3 result file,
#' although the order can be different.
#'
#' You can put any useful information in the sample information file.
#' Recommended columns are:
#' - `group`: grouping or conditions, e.g. "control" or "tumor",
#'   required for most downstream analyses
#' - `batch`: batch information, required for batch effect correction
#'
#' # Protein inference
#'
#' By default, this function automatically performs protein inference using the
#' parsimony method to resolve shared glycopeptides. This converts the plural
#' columns (`proteins`, `genes`, `protein_sites`) to singular equivalents
#' (`protein`, `gene`, `protein_site`).
#' 
#' # Aggregation
#'
#' pGlycoQuant performs quantification on the PSM level.
#' This level of information is too detailed for most downstream analyses.
#' This function aggregate PSMs into glycopeptides through summation.
#' For each glycopeptide (unique combination of "peptide", "peptide_site", "protein", "protein_site",
#' "gene", "glycan_composition", "glycan_structure"),
#' we sum up the quantifications of all PSMs that belong to this glycopeptide.
#'
#' # Output
#'
#' This function returns a [glyexp::experiment()] object.
#'
#' The following columns could be found in the variable information tibble:
#' - `peptide`: character, peptide sequence
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accession (after protein inference)
#' - `protein_site`: integer, site of glycosylation on protein (after protein inference)
#' - `gene`: character, gene name (symbol) (after protein inference)
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures.
#'
#' @param fp File path of the pGlycoQuant result file.
#' @param sample_info File path of the sample information file (csv),
#'  or a sample information data.frame/tibble.
#' @param quant_method Quantification method. Either "label-free" or "TMT".
#' @param glycan_type Glycan type. Either "N" or "O". Default is "N".
#' @param sample_name_converter A function to convert sample names from file paths.
#'  The function should take a character vector of old sample names
#'  and return new sample names.
#'  Note that sample names in `sample_info` should match the new names.
#'  If NULL, original names are kept.
#'
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()],
#'   [glyrepr::glycan_structure()]
#' @export
read_pglyco3_pglycoquant <- function(
  fp,
  sample_info = NULL,
  quant_method = c("label-free", "TMT"),
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".list")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)

  # ----- Read data -----
  # Keep all PSM-level columns for variable info
  if (quant_method == "label-free") {
    exp <- .read_pglyco3_pglycoquant_label_free(
      fp,
      sample_info,
      glycan_type,
      sample_name_converter
    )
  } else {
    rlang::abort("TMT quantification is not supported yet.")
  }

  exp
}


.read_pglyco3_pglycoquant_label_free <- function(
  fp,
  sample_info = NULL,
  glycan_type,
  sample_name_converter
) {
  # ----- Read data -----
  cli::cli_progress_step("Reading data")
  df <- .read_pglyco3_file_into_tibble(fp)
  expr_mat <- .extract_expr_mat_from_pglycoquant(df)
  samples <- .extract_sample_names_from_pglycoquant(df, sample_name_converter)
  sample_info <- .process_sample_info(sample_info, samples, glycan_type)
  var_info <- .extract_var_info_from_pglyco3(df)
  colnames(expr_mat) <- sample_info$sample
  rownames(expr_mat) <- var_info$variable

  # ----- Protein inference -----
  cli::cli_progress_step("Performing protein inference")
  var_info <- .infer_proteins(var_info)

  # ----- Aggregate PSMs to glycopeptides -----
  cli::cli_progress_step("Aggregating PSMs to glycopeptides")
  aggregated_result <- .aggregate_wide(var_info, expr_mat)
  var_info <- aggregated_result$var_info
  expr_mat <- aggregated_result$expr_mat

  # ----- Parse glycan compositions and structures -----
  cli::cli_progress_step("Parsing glycan compositions and structures")
  var_info <- dplyr::mutate(
    var_info,
    glycan_composition = .convert_pglyco3_comp(.data$glycan_composition),
    glycan_structure = glyparse::parse_pglyco_struc(.data$glycan_structure)
  )

  # ----- Pack Experiment -----
  cli::cli_progress_step("Packing experiment")
  glyexp::experiment(
    expr_mat, sample_info, var_info,
    exp_type = "glycoproteomics",
    glycan_type = glycan_type,
    quant_method = "label-free"
  )
}


# Extract variable information from pGlyco3 result
# Glycan composition and structure are not parsed here
# A "variable" column is added to the tibble
.extract_var_info_from_pglyco3 <- function(df) {
  new_names <- c(
    peptide = "Peptide",
    proteins = "Proteins",
    genes = "Genes",
    glycan_composition = "GlycanComposition",
    glycan_structure = "PlausibleStruct",
    peptide_site = "GlySite",
    protein_sites = "ProSites"
  )
  df %>%
    dplyr::select(all_of(new_names)) %>%
    dplyr::mutate(
      genes = stringr::str_remove(.data$genes, ";$"),
      proteins = stringr::str_replace_all(.data$proteins, "sp\\|(\\w+)\\|\\w+", "\\1")
    ) %>%
    dplyr::mutate(variable = stringr::str_c("PSM", dplyr::row_number()), .before = 1)
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

  # TODO: check column existence
  suppressWarnings(
    suppressMessages(readr::read_tsv(fp, col_types = col_types, progress = FALSE)),
    classes = "vroom_mismatched_column_name"
  )
}


.convert_pglyco3_comp <- function(x) {
  # Define mapping from pGlyco3 notation to generic monosaccharides
  pglyco_to_generic <- c(
    "H" = "Hex",      # Hexose
    "N" = "HexNAc",   # N-Acetylhexosamine  
    "A" = "NeuAc",    # N-Acetylneuraminic acid
    "G" = "HexA",     # Hexuronic acid
    "F" = "dHex"      # Deoxyhexose (Fucose)
  )
  
  extract_n_mono <- function(comp, mono) {
    n <- stringr::str_extract(comp, paste0(mono, "\\((\\d+)\\)"), group = 1)
    dplyr::if_else(is.na(n), 0L, as.integer(n))
  }

  unique_x <- unique(x)
  comp_df <- purrr::map_dfc(names(pglyco_to_generic), ~ {
    counts <- purrr::map_int(unique_x, extract_n_mono, mono = .x)
    tibble::tibble(!!pglyco_to_generic[[.x]] := counts)
  })
  
  # Convert each row to a glyrepr_composition object for unique values
  unique_compositions <- purrr::pmap(comp_df, function(...) {
    counts <- c(...)
    counts <- counts[counts > 0]
    if (length(counts) == 0) {
      return(glyrepr::as_glycan_composition(integer(0)))
    }
    glyrepr::as_glycan_composition(counts)
  })

  unique_compositions <- do.call(c, unique_compositions)
  compositions <- unique_compositions[match(x, unique_x)]
  compositions
}

# ----- Internal protein inference functions -----
.infer_proteins <- function(var_info) {
  # Parse the proteins, genes, and protein_sites columns
  proteins_list <- stringr::str_split(var_info$proteins, ";")
  genes_list <- stringr::str_split(var_info$genes, ";")
  protein_sites_list <- stringr::str_split(var_info$protein_sites, ";")

  # Create a mapping from protein to glycopeptide indices
  protein_to_gp <- proteins_list %>%
    purrr::imap(~ tibble::tibble(protein = .x, gp_idx = .y)) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(.data$protein) %>%
    dplyr::summarise(gp_indices = list(.data$gp_idx), .groups = "drop") %>%
    (function(x) stats::setNames(x$gp_indices, x$protein))
  
  # Greedy set cover algorithm
  uncovered_gps <- seq_len(nrow(var_info))
  selected_proteins <- character(0)
  
  while (length(uncovered_gps) > 0) {
    # Find the protein that covers the most uncovered glycopeptides
    available_proteins <- setdiff(names(protein_to_gp), selected_proteins)
    
    if (length(available_proteins) == 0) break
    
    # Calculate coverage for all available proteins
    coverages <- available_proteins %>%
      purrr::map_int(~ length(intersect(protein_to_gp[[.x]], uncovered_gps))) %>%
      purrr::set_names(available_proteins)
    
    # Find the protein with maximum coverage
    max_coverage <- max(coverages)
    if (max_coverage == 0) break
    
    best_protein <- names(coverages)[which.max(coverages)]
    
    # Add the best protein to the selected set
    selected_proteins <- c(selected_proteins, best_protein)
    
    # Remove covered glycopeptides from uncovered set
    uncovered_gps <- setdiff(uncovered_gps, protein_to_gp[[best_protein]])
  }
  
  # Calculate coverage count for each selected protein
  protein_coverage_count <- selected_proteins %>%
    purrr::map_int(~ length(protein_to_gp[[.x]])) %>%
    purrr::set_names(selected_proteins)
  
  # Assign each glycopeptide to the selected protein with most coverage
  assignments <- purrr::pmap(
    list(proteins_list, genes_list, protein_sites_list),
    function(prots, genes, sites) {
      # Find available selected proteins for this glycopeptide
      available_selected <- intersect(prots, selected_proteins)
      
      if (length(available_selected) > 0) {
        # Choose the protein with highest coverage count
        # If tie, choose the first one
        best_protein <- available_selected[which.max(protein_coverage_count[available_selected])]
        
        # Find the corresponding gene and site
        best_idx <- which(prots == best_protein)[1]
        
        list(
          protein = best_protein,
          gene = genes[best_idx],
          site = sites[best_idx]
        )
      } else {
        list(protein = NA_character_, gene = NA_character_, site = NA_character_)
      }
    }
  )
  
  assigned_protein <- purrr::map_chr(assignments, ~ .x$protein)
  assigned_gene <- purrr::map_chr(assignments, ~ .x$gene) 
  assigned_site <- purrr::map_chr(assignments, ~ .x$site)
  
  # Update var_info with protein inference results
  new_var_info <- var_info %>%
    dplyr::select(-tidyselect::all_of(c("proteins", "genes", "protein_sites"))) %>%
    dplyr::mutate(
      protein = assigned_protein,
      gene = assigned_gene,
      protein_site = as.integer(assigned_site)
    )

  new_var_info
}