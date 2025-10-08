# Process sample information for read functions
# This is a common utility function for all read_* functions
.process_sample_info <- function(sample_info, samples, glycan_type) {
  if (is.null(sample_info)) {
    sample_info <- tibble::tibble(sample = samples)
    cli::cli_alert_info(
      "No sample information file passed in. An empty tibble will be used."
    )
  } else {
    if (is.character(sample_info)) {
      # Read from file
      sample_info <- suppressMessages(readr::read_csv(sample_info))
    } else if (is.data.frame(sample_info)) {
      # Convert data.frame to tibble
      sample_info <- tibble::as_tibble(sample_info)
    }  # Otherwise, assume it's already a tibble.
    
    # Validate sample_info format by creating a dummy experiment
    # This passes the checking responsibility to `experiment()`.
    local({
      fake_var_info <- tibble::tibble(
        variable = character(0),
        protein = character(0),
        protein_site = integer(0),
        glycan_composition = glyrepr::glycan_composition(),
      )
      fake_expr_mat <- matrix(
        nrow = 0, ncol = length(samples),
        dimnames = list(NULL, samples)
      )

      # This line will throw an error if sample_info is not in correct format.
      glyexp::experiment(
        fake_expr_mat, sample_info, fake_var_info,
        exp_type = "glycoproteomics",
        glycan_type = glycan_type
      )
    })
  }
  
  sample_info
}
