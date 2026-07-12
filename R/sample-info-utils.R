#' Process sample information for read functions
#'
#' @param sample_info A file path, data frame, or `NULL`.
#' @param samples The sample identifiers present in the input data.
#'
#' @returns A validated sample information tibble.
#' @noRd
.process_sample_info <- function(sample_info, samples) {
  samples <- as.character(samples)

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
    } # Otherwise, assume it's already a tibble.
  }

  checkmate::assert_names(colnames(sample_info), must.include = "sample")
  sample_info$sample <- as.character(sample_info$sample)
  checkmate::assert_character(
    sample_info$sample,
    any.missing = FALSE,
    unique = TRUE
  )
  if (!setequal(sample_info$sample, samples)) {
    cli::cli_abort(
      "Sample identifiers in {.arg sample_info} must match the input data."
    )
  }

  factor_cols <- intersect(
    c("group", "batch", "bio_rep"),
    colnames(sample_info)
  )
  sample_info[factor_cols] <- purrr::map(sample_info[factor_cols], factor)

  sample_info
}
