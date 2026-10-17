test_that("protein groups survive inference in each reader", {
  df <- tibble::tibble(
    peptide = c("NAT", "NVT"),
    proteins = c("P2;P1", "P1"),
    genes = c("G2;G1", "G1"),
    protein_sites = c("2;1", "3")
  )
  inferred <- suppressMessages(.infer_proteins_df(df))
  expect_identical(inferred$protein_group, df$proteins)
  expect_identical(inferred$protein, c("P1", "P1"))

  sites <- .resolve_glyco_decipher_sites(c("P2@2;P1@1;P1@3", "P1@4"))
  expect_identical(sites$protein_group, c("P2;P1", "P1"))
  expect_identical(sites$protein, c("P1", "P1"))
  expect_identical(sites$protein_site, c(NA_integer_, 4L))

  strucgp <- tibble::tibble(
    file_name = "S1",
    peptide = c("NAT", "NVT"),
    protein_id = df$proteins,
    gene_name = df$genes,
    glycosite_position = df$protein_sites,
    glycan_composition = "H3N2",
    structure_coding = "unused"
  )
  tidy <- suppressMessages(.tidy_strucgp(strucgp, parse_structure = FALSE))
  expect_identical(tidy$protein_group, df$proteins)
  expect_identical(tidy$protein, c("P1", "P1"))
})

test_that("group annotations preserve aggregation and align across samples", {
  df <- tibble::tibble(
    peptide = "NAT",
    protein = "P1",
    protein_site = 1L,
    glycan_composition = "HexNAc(1)",
    protein_group = c("P2;P1", "P1;P3", "P1;P4"),
    sample = c("S1", "S1", "S2"),
    value = c(10, 20, 40)
  )
  read <- function(x) {
    suppressMessages(.read_template(
      x,
      NULL,
      "N",
      "label-free",
      composition_parser = glyrepr::as_glycan_composition,
      parse_structure = FALSE
    ))
  }
  result <- read(df)
  baseline <- read(dplyr::select(df, -"protein_group"))
  expect_identical(
    SummarizedExperiment::assay(result),
    SummarizedExperiment::assay(baseline)
  )
  expect_identical(
    SummarizedExperiment::rowData(result)$protein_group,
    "P2;P1;P3;P4"
  )
  expect_equal(
    unname(SummarizedExperiment::assay(result)),
    matrix(c(30, 40), 1)
  )
})

test_that("all protein inference readers retain groups in rowData", {
  readers <- list(
    function() {
      read_pglyco3_pglycoquant(test_path(
        "data/pglyco3-pglycoquant-LFQ-result.list"
      ))
    },
    function() read_pglyco3(test_path("data/pglyco3-LFQ-result.txt")),
    function() read_glyco_decipher(test_path("data/glyco-decipher-result.csv")),
    function() {
      read_strucgp(
        test_path("data/strucgp-result.xlsx"),
        parse_structure = FALSE
      )
    }
  )
  for (read in readers) {
    result <- suppressMessages(read())
    annotations <- SummarizedExperiment::rowData(result)
    expect_type(annotations$protein_group, "character")
    expect_length(annotations$protein_group, nrow(result))
    expect_equal(
      vapply(
        seq_len(nrow(result)),
        function(i) {
          annotations$protein[i] %in%
            strsplit(annotations$protein_group[i], ";", fixed = TRUE)[[1]]
        },
        logical(1)
      ),
      rep(TRUE, nrow(result))
    )
  }
})
