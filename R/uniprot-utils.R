# ----- UniProt identifier parsing utilities -----

# Internal function to extract UniProt accession from protein identifier strings
#
# This function handles various UniProt identifier formats including:
# - SwissProt format: "sp|P08185|CBG_HUMAN" -> "P08185"
# - TrEMBL format: "tr|A0A024R6P0|A0A024R6P0_HUMAN" -> "A0A024R6P0"
# - Isoform format: "sp|P08185-1|CBG_HUMAN" -> "P08185-1"
# - With prefix: ">sp|P08185|CBG_HUMAN" -> "P08185"
# - Multiple proteins: "sp|P08185|CBG_HUMAN;sp|P19652|A1AG2_HUMAN" -> "P08185;P19652"
# - Malformed entries return NA
#
# protein_ids: A character vector of protein identifier strings
# Returns: A character vector of extracted UniProt accessions
.extract_uniprot_accession <- function(protein_ids) {
  # First try to extract using str_replace_all for multiple proteins
  result <- stringr::str_replace_all(protein_ids, "(?:>)?(?:sp|tr)\\|([\\w-]+)\\|[^;]*", "\\1")

  # For entries that didn't match (result == original), assign NA
  no_match <- result == protein_ids
  result[no_match] <- NA_character_

  result
}

.add_gene_symbols <- function(df, orgdb) {
  if (!requireNamespace(orgdb, quietly = TRUE)) {
    cli::cli_alert_info("Package `{orgdb}` not installed. Skipping gene symbol conversion.")
    return(df)
  }

  # Get unique protein IDs (uniprot IDs)
  unique_proteins <- unique(df$protein)
  unique_proteins <- unique_proteins[!is.na(unique_proteins)]

  if (length(unique_proteins) == 0) {
    cli::cli_alert_info("No valid protein IDs found, skipping gene symbol conversion.")
    return(df)
  }

  # Try to convert uniprot IDs to gene symbols
  tryCatch(
    {
      org_db <- utils::getFromNamespace(orgdb, orgdb)

      # Use mapIds to convert UNIPROT to SYMBOL
      mapped_symbols <- suppressMessages(AnnotationDbi::mapIds(
        org_db,
        keys = unique_proteins,
        column = "SYMBOL",
        keytype = "UNIPROT",
        multiVals = "first"
      ))

      id_mapping <- tibble::tibble(
        UNIPROT = names(mapped_symbols),
        SYMBOL = mapped_symbols
      ) %>%
        tidyr::drop_na("SYMBOL")

      # Join with var_info to add gene symbols
      df <- dplyr::left_join(
        df,
        id_mapping,
        by = c("protein" = "UNIPROT")
      ) %>%
        dplyr::rename(gene = "SYMBOL")

      # Report conversion success
      n_converted <- sum(!is.na(unique(df$gene)))
      n_total <- length(unique(df$protein))
      perc_converted <- round(n_converted / n_total * 100, 1)
      cli::cli_alert_info(
        "Converted {.val {n_converted}} of {.val {n_total}} ({.val {perc_converted}}%) protein IDs to gene symbols."
      )
    },
    error = function(e) {
      cli::cli_alert_warning(
        "Failed to convert protein IDs to gene symbols: {.val {e$message}}"
      )
    }
  )

  df
}