# ----- Internal protein inference functions -----

#' Internal function for protein inference using parsimony method
#'
#' This function takes a character vector of proteins (with multiple proteins
#' separated by ";") and finds the minimal set of proteins that can explain
#' all observed entries using a greedy set cover algorithm.
#'
#' @param proteins_vec A character vector where each element contains one or more
#'   proteins separated by ";"
#' @returns An integer vector indicating the index of the selected protein for
#'   each entry in the input vector
#' @noRd
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

  # Map proteins to integer IDs for faster comparisons
  protein_names <- names(protein_to_entries)
  protein_entry_list <- protein_to_entries
  protein_entry_list <- lapply(protein_entry_list, unique)
  protein_entry_list <- unname(protein_entry_list)
  protein_id_lookup <- stats::setNames(seq_along(protein_names), protein_names)
  protein_rank <- integer(length(protein_names))
  protein_rank[order(protein_names)] <- seq_along(protein_names)

  # Greedy set cover algorithm
  entry_to_proteins <- lapply(proteins_list, function(prots) {
    ids <- protein_id_lookup[prots]
    ids <- ids[!is.na(ids)]
    unique(unname(ids))
  })
  uncovered <- rep_len(TRUE, length(proteins_vec))
  selected_protein_ids <- integer(0)
  cover <- lengths(protein_entry_list)
  storage.mode(cover) <- "double"
  
  while (any(uncovered)) {
    # Find the protein that covers the most uncovered entries
    max_coverage <- max(cover)
    if (!is.finite(max_coverage) || max_coverage <= 0) break
    
    # Find the protein with maximum coverage
    candidate_ids <- which(cover == max_coverage)
    if (length(candidate_ids) == 0) break
    best_protein_id <- candidate_ids[which.min(protein_rank[candidate_ids])]
    
    newly_covered <- protein_entry_list[[best_protein_id]]
    if (is.null(newly_covered) || length(newly_covered) == 0) {
      cover[best_protein_id] <- -Inf
      next
    }
    
    newly_covered <- newly_covered[uncovered[newly_covered]]
    if (length(newly_covered) == 0) {
      cover[best_protein_id] <- -Inf
      next
    }
    
    # Add the best protein to the selected set
    selected_protein_ids <- c(selected_protein_ids, best_protein_id)
    
    # Remove covered entries from uncovered set
    uncovered[newly_covered] <- FALSE
    
    # Decrement coverage counts for proteins associated with newly covered entries
    for (entry_idx in newly_covered) {
      entry_proteins <- entry_to_proteins[[entry_idx]]
      if (length(entry_proteins) == 0) next
      cover[entry_proteins] <- cover[entry_proteins] - 1
    }
    
    cover[cover < 0 & is.finite(cover)] <- 0
    cover[best_protein_id] <- -Inf
  }

  # Calculate coverage count for each selected protein
  selected_proteins <- protein_names[selected_protein_ids]
  protein_coverage_count <- selected_protein_ids %>%
    purrr::map_int(~ length(protein_entry_list[[.x]])) %>%
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

#' Perform protein inference on a list of vectors
#'
#' This function takes a list of vectors (proteins, genes, protein_sites) and
#' performs protein inference on the first vector (proteins), then applies the
#' same selection to all other vectors in the list.
#'
#' @param protein_vectors A named list containing vectors for proteins, genes,
#'   and protein_sites
#' @returns A named list with the same structure but with inferred selections
#' @noRd
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

#' Wrapper function for data frame input (for backward compatibility)
#'
#' This function maintains the original interface for pglyco3 functions
#'
#' @param df A data frame containing the columns "peptide", "proteins", "genes",
#'   and "protein_sites"
#' @returns A data frame with the same structure but with inferred selections
#' @noRd
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
