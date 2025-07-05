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
#' Three methods are available:
#' - `parsimony`: Uses a greedy set cover algorithm to find the minimal set of
#'   proteins that explain all glycopeptides, then assigns each shared glycopeptide
#'   to the protein with highest coverage.
#' - `unique`: Only retains glycopeptides that are uniquely assigned to a single
#'   protein.
#' - `share`: Splits shared glycopeptides into multiple entries, with expression
#'   values divided equally among associated proteins.
#'
#' # Output
#'
#' This function returns a [glyexp::experiment()] object.
#'
#' The following columns could be found in the variable information tibble:
#' - `charge`: integer, charge state
#' - `peptide`: character, peptide sequence
#' - `modifications`: character, modifications other than glycosylation,
#'   separated by semicolon, e.g. `5,Carbamidomethyl[C];10,Carbamidomethyl[C]`
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' - `glycan_structure`: [glyrepr::glycan_structure()], glycan structures.
#' - `peptide_site`: integer, site of glycosylation on peptide
#' - `protein`: character, protein accession (after protein inference)
#' - `gene`: character, gene name (symbol) (after protein inference)
#' - `protein_site`: integer, site of glycosylation on protein (after protein inference)
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
#' @param protein_inference_method Method for protein inference. Either "parsimony",
#'  "unique", or "share". Default is "parsimony".
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
  sample_name_converter = NULL,
  protein_inference_method = c("parsimony", "unique", "share")
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".list")
  .check_sample_info_arg(sample_info)
  quant_method <- rlang::arg_match(quant_method, c("label-free", "TMT"))
  .check_sample_name_conv_arg(sample_name_converter)
  glycan_type <- rlang::arg_match(glycan_type)
  protein_inference_method <- rlang::arg_match(protein_inference_method)

  # ----- Read data -----
  # Keep all PSM-level columns for variable info
  if (quant_method == "label-free") {
    exp <- .read_pglyco3_pglycoquant_label_free(
      fp,
      sample_info,
      glycan_type,
      sample_name_converter,
      protein_inference_method
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
  sample_name_converter,
  protein_inference_method
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

  # ----- Parse glycan compositions and structures -----
  cli::cli_progress_step("Parsing glycan compositions and structures")
  var_info <- dplyr::mutate(
    var_info,
    glycan_composition = .convert_pglyco3_comp(.data$glycan_composition),
    glycan_structure = glyparse::parse_pglyco_struc(.data$glycan_structure)
  )

  # ----- Pack Experiment -----
  cli::cli_progress_step("Packing experiment")
  exp <- glyexp::experiment(
    expr_mat, sample_info, var_info,
    exp_type = "glycoproteomics",
    glycan_type = glycan_type,
    quant_method = "label-free"
  )

  # ----- Protein inference -----
  cli::cli_progress_step("Performing protein inference")
  exp <- .infer_protein_internal(exp, protein_inference_method)

  exp
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
    protein_sites = "ProSites",
    charge = "Charge",
    modifications = "Mod"
  )
  df %>%
    dplyr::select(all_of(new_names)) %>%
    dplyr::mutate(
      modifications = stringr::str_remove(.data$modifications, ";$"),
      modifications = dplyr::if_else(is.na(.data$modifications), "", .data$modifications),
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

.infer_protein_internal <- function(exp, method) {
  if (!"proteins" %in% colnames(exp$var_info)) {
    return(exp)  # No protein inference needed
  }
  
  switch(method,
    unique = .infer_pro_unique_internal(exp),
    parsimony = .infer_pro_parsimony_internal(exp),
    share = .infer_pro_share_internal(exp)
  )
}

.infer_pro_unique_internal <- function(exp) {
  exp %>%
    glyexp::filter_var(!stringr::str_detect(.data$proteins, ";")) %>%
    glyexp::rename_var(tidyselect::all_of(c("protein" = "proteins", "gene" = "genes", "protein_site" = "protein_sites"))) %>%
    glyexp::mutate_var(protein_site = as.integer(.data$protein_site))
}

.infer_pro_parsimony_internal <- function(exp) {
  var_info <- exp$var_info
  
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
  
  # Update the experiment
  new_var_info <- var_info %>%
    dplyr::select(-tidyselect::all_of(c("proteins", "genes", "protein_sites"))) %>%
    dplyr::mutate(
      protein = assigned_protein,
      gene = assigned_gene,
      protein_site = as.integer(assigned_site)
    )
  
  new_exp <- exp
  new_exp$var_info <- new_var_info
  new_exp
}

.infer_pro_share_internal <- function(exp) {
  var_info <- exp$var_info
  expr_mat <- exp$expr_mat
  
  # Parse the proteins, genes, and protein_sites columns
  proteins_list <- stringr::str_split(var_info$proteins, ";")
  genes_list <- stringr::str_split(var_info$genes, ";")
  protein_sites_list <- stringr::str_split(var_info$protein_sites, ";")
  
  # Helper function to create a single protein entry from shared glycopeptide
  create_protein_entry <- function(protein, gene, site, idx, shared_expr) {
    new_variable <- paste0(var_info$variable[idx], "_", protein)
    
    # Create new var_info row
    new_var_row <- var_info[idx, , drop = FALSE]
    new_var_row$variable <- new_variable
    new_var_row <- new_var_row %>%
      dplyr::select(-tidyselect::all_of(c("proteins", "genes", "protein_sites"))) %>%
      dplyr::mutate(
        protein = protein,
        gene = gene,
        protein_site = as.integer(site)
      )
    
    # Create new expression row
    new_expr_row <- shared_expr
    rownames(new_expr_row) <- new_variable
    
    list(var_info = new_var_row, expr_row = new_expr_row)
  }
  
  # Helper function to process a single glycopeptide and split by proteins
  process_glycopeptide <- function(proteins, genes, sites, idx) {
    n_proteins <- length(proteins)
    original_expr <- expr_mat[idx, , drop = FALSE]
    shared_expr <- original_expr / n_proteins
    
    # Create rows for each protein
    purrr::pmap(
      list(protein = proteins, gene = genes, site = sites),
      ~ create_protein_entry(..1, ..2, ..3, idx, shared_expr)
    )
  }
  
  # Process each glycopeptide and split by proteins
  share_results <- purrr::pmap(
    list(
      proteins = proteins_list,
      genes = genes_list, 
      sites = protein_sites_list,
      idx = seq_len(nrow(var_info))
    ),
    process_glycopeptide
  ) %>%
    purrr::flatten()
  
  # Extract var_info and expression data
  new_var_info_list <- purrr::map(share_results, ~ .x$var_info)
  new_expr_rows <- purrr::map(share_results, ~ .x$expr_row)
  
  # Combine all new data
  new_var_info <- dplyr::bind_rows(new_var_info_list)
  new_expr_mat <- do.call(rbind, new_expr_rows)
  
  # Create new experiment
  new_exp <- exp
  new_exp$var_info <- new_var_info
  new_exp$expr_mat <- new_expr_mat
  new_exp
}
