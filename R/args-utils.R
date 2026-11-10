#' Validate common arguments for read_* functions
#'
#' This function provides centralized validation for parameters commonly used
#' across all read_* functions in glyread package. It reduces code duplication
#' and ensures consistent validation logic.
#'
#' @param sample_info Sample information parameter
#' @param quant_method Quantification method parameter
#' @param glycan_type Glycan type parameter
#' @param sample_name_converter Sample name converter function parameter
#' @param parse_structure Parse structure parameter (optional, only for some functions)
#' @param orgdb OrgDb parameter (optional, only for some functions)
#'
#' @noRd
.validate_read_args <- function(
  sample_info = NULL,
  quant_method = "label-free",
  glycan_type = "N",
  sample_name_converter = NULL,
  parse_structure = NULL,
  orgdb = NULL
) {
  if (!.is_http_url(sample_info)) {
    checkmate::assert(
      checkmate::check_null(sample_info),
      checkmate::check_file_exists(
        sample_info,
        access = "r",
        extension = ".csv"
      ),
      checkmate::check_data_frame(sample_info)
    )
  }

  # Validate and match quant_method
  checkmate::assert_choice(quant_method, c("label-free", "TMT"))

  # Validate and match glycan_type
  checkmate::assert_choice(
    glycan_type,
    c("N", "O-GalNAc", "O-GlcNAc", "O-Man", "O-Fuc", "O-Glc")
  )

  # Validate sample_name_converter
  checkmate::assert(
    checkmate::check_null(sample_name_converter),
    checkmate::check_function(sample_name_converter)
  )

  # Validate parse_structure if provided
  if (!is.null(parse_structure)) {
    checkmate::assert_logical(parse_structure, len = 1)
  }

  # Validate orgdb if provided
  if (!is.null(orgdb)) {
    checkmate::assert_string(orgdb)
  }
}

.is_http_url <- function(x) {
  is.character(x) &&
    length(x) == 1L &&
    !is.na(x) &&
    grepl("^https?://", x, ignore.case = TRUE)
}

.raw_file_url <- function(url) {
  sub(
    "^https?://github[.]com/([^/?#]+)/([^/?#]+)/blob/",
    "https://raw.githubusercontent.com/\\1/\\2/",
    url,
    ignore.case = TRUE
  )
}

.read_file_input <- function(fp, extensions) {
  if (.is_http_url(fp)) {
    return(.raw_file_url(fp))
  }
  checkmate::assert_file_exists(fp, access = "r", extension = extensions)
  fp
}
