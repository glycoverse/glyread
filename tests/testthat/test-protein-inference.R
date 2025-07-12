test_that(".infer_proteins_internal works with no ambiguity", {
  proteins <- c("PROA", "PROB", "PROC")
  result <- .infer_proteins_internal(proteins)
  expect_equal(result, c(1, 1, 1))
})

test_that(".infer_proteins_internal resolves basic ambiguity correctly", {
  proteins <- c("PROA", "PROB;PROA", "PROC")
  result <- .infer_proteins_internal(proteins)
  expect_equal(result, c(1, 2, 1))  # Should select PROA for second entry
})

test_that(".infer_proteins_internal handles proteins appearing in multiple entries", {
  proteins <- c("PROA;PROB", "PROA", "PROB")
  result <- .infer_proteins_internal(proteins)
  # PROA covers entries 1,2 (2 entries), PROB covers entries 1,3 (2 entries)
  # Algorithm should pick one of them first, then the other
  expect_length(result, 3)
  expect_true(all(result %in% c(1, 2)))
  expect_true(all(!is.na(result)))
})

test_that(".infer_proteins_internal handles complex parsimony cases", {
  # Test case: Parsimony principle - choose minimal set
  proteins <- c(
    "PROA;PROB;PROC",  # Entry 1: all three proteins
    "PROA;PROB",       # Entry 2: PROA and PROB
    "PROA",            # Entry 3: only PROA
    "PROB;PROC",       # Entry 4: PROB and PROC
    "PROC",            # Entry 5: only PROC
    "PROC"             # Entry 6: only PROC
  )
  result <- .infer_proteins_internal(proteins)

  # Firstly, the minimum set that covers all is PROA + PROC,
  # so the last four index should be 1, 1, 2, 1.
  # As for the first one, both PROA and PROC appear.
  # We should pick the one appears the most.
  # PROA appears 3 times, and PROC appears 4 times, so we should pick PROC.
  expected <- c(3, 1, 1, 2, 1, 1)
  expect_equal(result, expected)
})

test_that(".infer_proteins_internal handles empty input", {
  expect_equal(.infer_proteins_internal(character(0)), integer(0))
})

test_that(".infer_proteins_internal handles single entry", {
  result <- .infer_proteins_internal("PROA")
  expect_equal(result, 1)
})

test_that(".infer_proteins_internal handles single entry with multiple proteins", {
  result <- .infer_proteins_internal("PROA;PROB;PROC")
  expect_equal(result, 1)  # Should select first protein
})

test_that(".infer_proteins_internal handles identical proteins across entries", {
  proteins <- c("PROA", "PROA", "PROA")
  result <- .infer_proteins_internal(proteins)
  expect_equal(result, c(1, 1, 1))
})

test_that(".infer_proteins_internal handles tie-breaking consistently", {
  # When multiple proteins have the same coverage, should pick consistently
  proteins <- c("PROA;PROB", "PROA;PROC", "PROB;PROC")
  result <- .infer_proteins_internal(proteins)

  # PROA + PROB, PROA + PROC, PROB + PROC all satisfy the demand,
  # and PROA, PROB, PROC appears the same amount of times.
  # In this case, use the one that appears first in alphabetical order,
  # which is PROA + PROB.
  # For the first one, both PROA and PROB appear, we use the alphabetical order again,
  # so we should pick PROA.
  expected <- c(1, 1, 1)
  expect_equal(result, expected)
})

test_that(".infer_proteins produces deterministic results regardless of protein order", {
  # Test all 8 possible combinations of protein order in 3 entries
  # Each entry has 2 proteins, so 2^3 = 8 combinations
  # Due to alphabetical tie-breaking, PROA should always be selected when available
  test_cases <- list(
    # Case 1: PROA;PROB, PROA;PROC, PROB;PROC
    list(
      proteins = c("PROA;PROB", "PROA;PROC", "PROB;PROC"),
      genes = c("GENEA;GENEB", "GENEA;GENEC", "GENEB;GENEC"),
      expected_proteins = c("PROA", "PROA", "PROB")
    ),
    # Case 2: PROB;PROA, PROA;PROC, PROB;PROC
    list(
      proteins = c("PROB;PROA", "PROA;PROC", "PROB;PROC"),
      genes = c("GENEB;GENEA", "GENEA;GENEC", "GENEB;GENEC"),
      expected_proteins = c("PROA", "PROA", "PROB")
    ),
    # Case 3: PROA;PROB, PROC;PROA, PROB;PROC
    list(
      proteins = c("PROA;PROB", "PROC;PROA", "PROB;PROC"),
      genes = c("GENEA;GENEB", "GENEC;GENEA", "GENEB;GENEC"),
      expected_proteins = c("PROA", "PROA", "PROB")
    ),
    # Case 4: PROB;PROA, PROC;PROA, PROB;PROC
    list(
      proteins = c("PROB;PROA", "PROC;PROA", "PROB;PROC"),
      genes = c("GENEB;GENEA", "GENEC;GENEA", "GENEB;GENEC"),
      expected_proteins = c("PROA", "PROA", "PROB")
    ),
    # Case 5: PROA;PROB, PROA;PROC, PROC;PROB
    list(
      proteins = c("PROA;PROB", "PROA;PROC", "PROC;PROB"),
      genes = c("GENEA;GENEB", "GENEA;GENEC", "GENEC;GENEB"),
      expected_proteins = c("PROA", "PROA", "PROB")
    ),
    # Case 6: PROB;PROA, PROA;PROC, PROC;PROB
    list(
      proteins = c("PROB;PROA", "PROA;PROC", "PROC;PROB"),
      genes = c("GENEB;GENEA", "GENEA;GENEC", "GENEC;GENEB"),
      expected_proteins = c("PROA", "PROA", "PROB")
    ),
    # Case 7: PROA;PROB, PROC;PROA, PROC;PROB
    list(
      proteins = c("PROA;PROB", "PROC;PROA", "PROC;PROB"),
      genes = c("GENEA;GENEB", "GENEC;GENEA", "GENEC;GENEB"),
      expected_proteins = c("PROA", "PROA", "PROB")
    ),
    # Case 8: PROB;PROA, PROC;PROA, PROC;PROB
    list(
      proteins = c("PROB;PROA", "PROC;PROA", "PROC;PROB"),
      genes = c("GENEB;GENEA", "GENEC;GENEA", "GENEC;GENEB"),
      expected_proteins = c("PROA", "PROA", "PROB")
    )
  )

  # Test each case
  for (i in seq_along(test_cases)) {
    case <- test_cases[[i]]
    protein_vectors <- list(
      proteins = case$proteins,
      genes = case$genes
    )

    result <- .infer_proteins(protein_vectors)

    # The selected proteins should always be the same regardless of order
    # PROA should always be selected due to alphabetical tie-breaking
    expect_equal(result$proteins, case$expected_proteins,
                 info = paste("Failed for case", i, "with proteins:", paste(case$proteins, collapse = ", ")))
  }
})

test_that(".infer_proteins works with basic list input", {
  protein_vectors <- list(
    proteins = c("PROA", "PROB;PROA", "PROC"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC"),
    protein_sites = c("123", "456;123", "789")
  )

  result <- .infer_proteins(protein_vectors)

  expected <- list(
    proteins = c("PROA", "PROA", "PROC"),
    genes = c("GENEA", "GENEA", "GENEC"),
    protein_sites = c("123", "123", "789")
  )
  expect_equal(result, expected)
})

test_that(".infer_proteins returns correct list structure", {
  protein_vectors <- list(
    proteins = c("PROA", "PROB;PROA", "PROC"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC"),
    protein_sites = c("123", "456;123", "789")
  )

  result <- .infer_proteins(protein_vectors)

  expect_type(result, "list")
  expect_equal(names(result), c("proteins", "genes", "protein_sites"))
  expect_equal(length(result$proteins), 3)
  expect_equal(length(result$genes), 3)
  expect_equal(length(result$protein_sites), 3)
})

test_that(".infer_proteins handles mismatched vector lengths", {
  # Test case: Different length vectors (should not work)
  protein_vectors <- list(
    proteins = c("PROA", "PROB;PROA"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC")
  )
  expect_error(.infer_proteins(protein_vectors), "must have the same length")
})

test_that(".infer_proteins handles empty list", {
  result <- .infer_proteins(list())
  expect_equal(result, list())
})

test_that(".infer_proteins handles single vector", {
  protein_vectors <- list(proteins = c("PROA", "PROB;PROA"))
  result <- .infer_proteins(protein_vectors)
  expect_equal(result, list(proteins = c("PROA", "PROA")))
})

test_that(".infer_proteins handles empty vectors", {
  protein_vectors <- list(
    proteins = character(0),
    genes = character(0)
  )
  result <- .infer_proteins(protein_vectors)
  expect_equal(result$proteins, character(0))
  expect_equal(result$genes, character(0))
})

test_that(".infer_proteins_df removes original columns", {
  test_df <- data.frame(
    peptide = c("PEP1", "PEP2", "PEP3"),
    proteins = c("PROA", "PROB;PROA", "PROC"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC"),
    protein_sites = c("123", "456;123", "789"),
    stringsAsFactors = FALSE
  )

  suppressMessages(result_df <- .infer_proteins_df(test_df))

  expect_false("proteins" %in% colnames(result_df))
  expect_false("genes" %in% colnames(result_df))
  expect_false("protein_sites" %in% colnames(result_df))
})

test_that(".infer_proteins_df adds new columns", {
  test_df <- data.frame(
    peptide = c("PEP1", "PEP2", "PEP3"),
    proteins = c("PROA", "PROB;PROA", "PROC"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC"),
    protein_sites = c("123", "456;123", "789"),
    stringsAsFactors = FALSE
  )

  suppressMessages(result_df <- .infer_proteins_df(test_df))

  expect_true("protein" %in% colnames(result_df))
  expect_true("gene" %in% colnames(result_df))
  expect_true("protein_site" %in% colnames(result_df))
})

test_that(".infer_proteins_df preserves other columns", {
  test_df <- data.frame(
    peptide = c("PEP1", "PEP2", "PEP3"),
    proteins = c("PROA", "PROB;PROA", "PROC"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC"),
    protein_sites = c("123", "456;123", "789"),
    other_col = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )

  suppressMessages(result_df <- .infer_proteins_df(test_df))

  expect_true("peptide" %in% colnames(result_df))
  expect_true("other_col" %in% colnames(result_df))
})

test_that(".infer_proteins_df produces correct data types", {
  test_df <- data.frame(
    peptide = c("PEP1", "PEP2", "PEP3"),
    proteins = c("PROA", "PROB;PROA", "PROC"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC"),
    protein_sites = c("123", "456;123", "789"),
    stringsAsFactors = FALSE
  )

  suppressMessages(result_df <- .infer_proteins_df(test_df))

  expect_type(result_df$protein_site, "integer")
})

test_that(".infer_proteins_df produces correct values", {
  test_df <- data.frame(
    peptide = c("PEP1", "PEP2", "PEP3"),
    proteins = c("PROA", "PROB;PROA", "PROC"),
    genes = c("GENEA", "GENEB;GENEA", "GENEC"),
    protein_sites = c("123", "456;123", "789"),
    stringsAsFactors = FALSE
  )

  suppressMessages(result_df <- .infer_proteins_df(test_df))

  expect_equal(result_df$protein, c("PROA", "PROA", "PROC"))
  expect_equal(result_df$gene, c("GENEA", "GENEA", "GENEC"))
  expect_equal(result_df$protein_site, c(123L, 123L, 789L))
})
