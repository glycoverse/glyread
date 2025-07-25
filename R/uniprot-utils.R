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
