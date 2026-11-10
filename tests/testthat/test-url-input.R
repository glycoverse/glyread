test_that("GitHub file pages resolve to raw data and direct URLs pass through", {
  expect_identical(
    .read_file_input(
      "https://github.com/org/repo/blob/main/data/a%20b.list",
      ".list"
    ),
    "https://raw.githubusercontent.com/org/repo/main/data/a%20b.list"
  )
  direct <- "https://example.org/org/repo/blob/main/download?token=abc"
  expect_identical(.read_file_input(direct, ".list"), direct)
})

test_that("HTTP results and sample information match their local files", {
  skip_if_not_installed("httpuv")
  root <- normalizePath(test_path("data"))
  port <- tryCatch(
    httpuv::randomPort(),
    httpuv_unavailable_port = function(e) skip("No local port is available")
  )
  server <- httpuv::startServer(
    "127.0.0.1",
    port,
    list(staticPaths = list("/" = root))
  )
  on.exit(server$stop(), add = TRUE)
  url <- function(file) paste0("http://127.0.0.1:", port, "/", file)

  result <- "pglyco3-pglycoquant-LFQ-result.list"
  samples <- "pglyco3-LFQ-sample-info.csv"
  expect_equal(
    read_pglyco3_pglycoquant(url(result), url(samples)),
    read_pglyco3_pglycoquant(file.path(root, result), file.path(root, samples))
  )

  workbook <- "strucgp-result.xlsx"
  expect_equal(
    read_strucgp(url(workbook), parse_structure = FALSE),
    read_strucgp(file.path(root, workbook), parse_structure = FALSE)
  )

  files <- file.path("msfragger_result", c("H_1", "H_2", "H_3"), "psm.tsv")
  expect_equal(
    read_msfragger(paste0(url(files), "?source=a/b")),
    read_msfragger(file.path(root, "msfragger_result"))
  )
})
