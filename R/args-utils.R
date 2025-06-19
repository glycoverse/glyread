.check_sample_info_arg <- function(sample_info) {
  checkmate::assert(
    checkmate::check_null(sample_info),
    checkmate::check_file_exists(
      sample_info, access = "r", extension = ".csv"
    ),
    checkmate::check_data_frame(sample_info)
  )
}

.check_sample_name_conv_arg <- function(sample_name_converter) {
  checkmate::assert(
    checkmate::check_null(sample_name_converter),
    checkmate::check_function(sample_name_converter)
  )
}
