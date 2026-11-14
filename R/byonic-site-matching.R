# Monoisotopic increments on the peptide, without reducing-end water.
.byonic_residue_mass <- c(
  Hex = 162.05282342,
  HexNAc = 203.07937252,
  dHex = 146.05790880,
  Pen = 132.04225873,
  HexA = 176.03208798,
  HexN = 161.06880784,
  dHexNAc = 187.08445790,
  ddHex = 130.06299418,
  Hep = 192.06338810,
  NeuAc = 291.09541651,
  NeuGc = 307.09033113,
  gKdn = 250.06886741,
  P = 79.96633089,
  S = 79.95681503,
  Ac = 42.01056468,
  Me = 14.01565006
)

.byonic_glycan_lookup <- function(source) {
  parsed <- lapply(source, function(x) {
    tryCatch(glyrepr::as_glycan_composition(x), error = function(e) NULL)
  })
  composition <- vapply(
    parsed,
    function(x) {
      if (is.null(x) || is.na(x)) NA_character_ else as.character(x)
    },
    character(1)
  )
  mass <- vapply(
    parsed,
    function(x) {
      if (is.null(x) || is.na(x)) {
        return(NA_real_)
      }
      counts <- as.list(glyrepr::convert_to_generic(x))[[1L]]
      sum(counts * .byonic_residue_mass[names(counts)])
    },
    numeric(1)
  )
  data.frame(source, composition, mass)
}

.byonic_site_pairs <- function(composition = character(), site = integer()) {
  data.frame(glycan_composition = composition, peptide_site = site)
}

.pair_byonic_glycans <- function(peptide, glycans, masses) {
  # Count residues, keeping all successive bracketed modifications at the
  # preceding site. An annotation before the first residue has no peptide site.
  core <- .extract_byonic_pglycoquant_peptide_core(peptide)
  tokens <- stringr::str_extract_all(core, "\\[[^]]*\\]|[A-Za-z]")[[1L]]
  tokens <- tokens[!is.na(tokens)]
  bracket <- startsWith(tokens, "[")
  positions <- cumsum(!bracket)
  annotated <- bracket & positions > 0L
  site <- positions[annotated]
  shifts <- suppressWarnings(as.numeric(stringr::str_sub(
    tokens[annotated],
    2L,
    -2L
  )))
  candidates <- lapply(masses, function(mass) {
    unique(site[which(abs(shifts - mass) <= 0.001)])
  })

  # Common path: each glycan has at most one candidate and no site is shared.
  hits <- as.integer(unlist(candidates))
  if (all(lengths(candidates) <= 1L) && !anyDuplicated(hits)) {
    sites <- vapply(candidates, \(x) x[1L], integer(1))
    return(unique(.byonic_site_pairs(glycans, sites)))
  }

  sites <- sort(unique(hits))
  edges <- do.call(rbind, lapply(candidates, \(x) sites %in% x))
  .resolve_byonic_site_matches(glycans, sites, edges)
}

# Rows are individual reported glycans and columns are distinct peptide sites.
.byonic_matching <- function(edges, weighted = NULL) {
  graph <- igraph::graph_from_biadjacency_matrix(edges, weighted = weighted)
  igraph::max_bipartite_match(graph)
}

.resolve_byonic_site_matches <- function(glycans, sites, edges) {
  best <- .byonic_matching(edges)
  size <- best$matching_size
  assigned <- best$matching[nrow(edges) + seq_along(sites)]
  # A site is certain if removing it reduces the maximum number of matches.
  occupied <- which(!is.na(assigned))
  certain <- Filter(
    function(j) {
      without <- edges
      without[, j] <- FALSE
      .byonic_matching(without)$matching_size < size
    },
    occupied
  )

  composition <- vapply(
    certain,
    function(j) {
      glycan <- glycans[assigned[j]]
      without <- edges
      without[which(glycans == glycan), j] <- FALSE
      if (.byonic_matching(without)$matching_size < size) {
        glycan
      } else {
        NA_character_
      }
    },
    character(1)
  )
  known_sites <- .byonic_site_pairs(composition, sites[certain])
  n_unknown <- length(glycans) - length(certain)
  if (n_unknown == 0L) {
    return(known_sites)
  }

  # At unknown sites retain each composition as often as it must occur there.
  # Among maximum matches, maximize its copies placed on the certain sites.
  remaining <- lapply(unique(glycans), function(glycan) {
    copies <- which(glycans %in% glycan)
    placed <- 0L
    if (any(edges[copies, certain, drop = FALSE])) {
      base_weight <- size + 1L
      weights <- edges * base_weight
      weights[copies, certain] <- weights[copies, certain] +
        edges[copies, certain]
      placed <- .byonic_matching(weights, weighted = TRUE)$matching_weight -
        size * base_weight
    }
    rep(glycan, length(copies) - placed)
  })
  remaining <- as.character(unlist(remaining))
  unknown <- c(remaining, rep(NA_character_, n_unknown - length(remaining)))
  unique(rbind(
    known_sites,
    .byonic_site_pairs(unknown, rep(NA_integer_, n_unknown))
  ))
}
