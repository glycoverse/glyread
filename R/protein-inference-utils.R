# ----- Internal protein inference functions -----

# Internal function for protein inference using parsimony method
#
# This function takes a character vector of proteins (with multiple proteins
# separated by ";") and finds the minimal set of proteins that can explain
# all observed entries using a greedy set cover algorithm.
#
# proteins_vec: A character vector where each element contains one or more
#   proteins separated by ";"
# Returns: An integer vector indicating the index of the selected protein for
#   each entry in the input vector
.infer_proteins_internal <- function(proteins_vec) {
  # Handle empty input
  if (length(proteins_vec) == 0) {
    return(integer(0))
  }

  # Parse the proteins
  proteins_list <- stringr::str_split(proteins_vec, ";")

  # Create a mapping from protein to entry indices
  protein_to_entries <- proteins_list %>%
    purrr::imap(~ tibble::tibble(protein = .x, entry_idx = .y)) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(.data$protein) %>%
    dplyr::summarise(entry_indices = list(.data$entry_idx), .groups = "drop") %>%
    (function(x) stats::setNames(x$entry_indices, x$protein))
  
  # Greedy set cover algorithm
  uncovered_entries <- seq_len(length(proteins_vec))
  selected_proteins <- character(0)
  
  while (length(uncovered_entries) > 0) {
    # Find the protein that covers the most uncovered entries
    available_proteins <- setdiff(names(protein_to_entries), selected_proteins)
    
    if (length(available_proteins) == 0) break
    
    # Calculate coverage for all available proteins
    coverages <- available_proteins %>%
      purrr::map_int(~ length(intersect(protein_to_entries[[.x]], uncovered_entries))) %>%
      purrr::set_names(available_proteins)
    
    # Find the protein with maximum coverage
    max_coverage <- max(coverages)
    if (max_coverage == 0) break

    # If tie, choose the first one in alphabetical order
    best_proteins <- names(coverages)[coverages == max_coverage]
    best_protein <- sort(best_proteins)[1]
    
    # Add the best protein to the selected set
    selected_proteins <- c(selected_proteins, best_protein)
    
    # Remove covered entries from uncovered set
    uncovered_entries <- setdiff(uncovered_entries, protein_to_entries[[best_protein]])
  }
  
  # Calculate coverage count for each selected protein
  protein_coverage_count <- selected_proteins %>%
    purrr::map_int(~ length(protein_to_entries[[.x]])) %>%
    purrr::set_names(selected_proteins)
  
  # Assign each entry to the selected protein with most coverage
  assignments <- purrr::map_int(proteins_list, function(prots) {
    # Find available selected proteins for this entry
    available_selected <- intersect(prots, selected_proteins)
    
    if (length(available_selected) > 0) {
      # Choose the protein with highest coverage count
      # If tie, choose the first one in alphabetical order
      max_count <- max(protein_coverage_count[available_selected])
      best_candidates <- available_selected[protein_coverage_count[available_selected] == max_count]
      best_protein <- sort(best_candidates)[1]
      
      # Find the index of the best protein in this entry's protein list
      which(prots == best_protein)[1]
    } else {
      # If no selected protein is available, return NA
      NA_integer_
    }
  })
  
  assignments
}

# Perform protein inference on a list of vectors
#
# This function takes a list of vectors (proteins, genes, protein_sites) and
# performs protein inference on the first vector (proteins), then applies the
# same selection to all other vectors in the list.
#
# protein_vectors: A named list containing vectors for proteins, genes,
#   and protein_sites
# Returns: A named list with the same structure but with inferred selections
.infer_proteins <- function(protein_vectors) {
  # Handle empty list
  if (length(protein_vectors) == 0) {
    return(list())
  }

  # Handle mismatched vector lengths
  if (dplyr::n_distinct(purrr::map_int(protein_vectors, length)) > 1) {
    cli::cli_abort("All vectors in `protein_vectors` must have the same length.")
  }

  # Perform protein inference on the proteins vector (first vector)
  first_vector_name <- names(protein_vectors)[1]
  selected_indices <- .infer_proteins_internal(protein_vectors[[first_vector_name]])
  
  # Apply the selection to all vectors in the list
  result <- purrr::map(protein_vectors, function(vec) {
    # Split each entry by ";"
    vec_list <- stringr::str_split(vec, ";")
    
    # Extract the selected element for each entry
    purrr::map2_chr(vec_list, selected_indices, function(elements, idx) {
      if (is.na(idx) || idx > length(elements)) NA_character_ else elements[idx]
    })
  })
  
  result
}

# Wrapper function for data frame input (for backward compatibility)
# This function maintains the original interface for pglyco3 functions
.infer_proteins_df <- function(df) {
  cli::cli_progress_step("Finding leader proteins")
  # Parse the proteins, genes, and protein_sites columns
  pep_df <- dplyr::distinct(df, .data$peptide, .data$proteins, .data$genes, .data$protein_sites)
  
  # Create the list of vectors for protein inference
  protein_vectors <- list(
    proteins = pep_df$proteins,
    genes = pep_df$genes,
    protein_sites = pep_df$protein_sites
  )
  
  # Perform protein inference using the new interface
  inferred_vectors <- .infer_proteins(protein_vectors)
  
  # Update var_info with protein inference results
  new_old_map <- pep_df %>%
    dplyr::mutate(
      protein = inferred_vectors$proteins,
      gene = inferred_vectors$genes,
      protein_site = as.integer(inferred_vectors$protein_sites)
    )
  new_df <- df %>%
    dplyr::left_join(new_old_map, by = c("peptide", "proteins", "genes", "protein_sites")) %>%
    dplyr::select(-tidyselect::all_of(c("proteins", "genes", "protein_sites")))

  new_df
}
