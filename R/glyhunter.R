#' Read GlyHunter result
#'
#' This function reads in a [GlyHunter](https://github.com/fubin1999/glyhunter)
#' result file and returns a [glyexp::experiment()] object.
#'
#' @details
#' # Which file to use?
#'
#' Use the "summary_area.csv" file in the GlyHunter result folder.
#' No edit is necessary.
#' 
#' # Variable information
#'
#' The following columns could be found in the variable information tibble:
#' - `glycan_composition`: [glyrepr::glycan_composition()], glycan compositions.
#' 
#' @inheritSection read_pglyco3_pglycoquant Sample information
#' 
#' @inheritParams read_pglyco3_pglycoquant
#' @param fp File path of the GlyHunter result file.
#' 
#' @returns An [glyexp::experiment()] object.
#' @seealso [glyexp::experiment()], [glyrepr::glycan_composition()]
#' @export
read_glyhunter <- function(
  fp,
  sample_info = NULL,
  glycan_type = c("N", "O"),
  sample_name_converter = NULL
) {
  # ----- Check arguments -----
  checkmate::assert_file_exists(fp, access = "r", extension = ".csv")
  .check_sample_info_arg(sample_info)
  glycan_type <- rlang::arg_match(glycan_type)
  .check_sample_name_conv_arg(sample_name_converter)

  # ----- Read data -----
  df <- suppressMessages(readr::read_csv(fp))

  # Prepare sample info
  samples <- setdiff(colnames(df), "glycan")
  sample_info <- .process_sample_info(sample_info, samples, glycan_type)

  # Prepare variable info
  var_info <- df %>%
    dplyr::select(all_of(c("variable" = "glycan"))) %>%
    dplyr::mutate(
      glycan_composition = stringr::str_remove_all(.data$variable, "\\[.*?\\]"),
      glycan_composition = glyrepr::as_glycan_composition(.data$glycan_composition)
    )

  # Prepare expression matrix
  expr_mat <- df %>%
    tibble::column_to_rownames("glycan") %>%
    as.matrix()

  # Pack experiment
  glyexp::experiment(
    expr_mat, sample_info, var_info,
    exp_type = "glycomics",
    glycan_type = glycan_type
  )
}
