.packed_backend <- function() {
  gsim:::.gsim_gbits_backend(Sys.getenv("GSIM_GBITS_LIBRARY"))
}

.packed_founder_fixture <- function(
  backend, n = 12L, seed = 2026, individual_offset = 0L,
  chromosome = "chrZ", return_genotypes = TRUE, return_segments = TRUE
) {
  marker_count <- 17L
  marker <- seq_len(marker_count)
  h1 <- outer(seq_len(4L), marker, function(donor, at) {
    (donor * 3L + at * 5L + at %/% 3L) %% 2L
  })
  h2 <- outer(seq_len(4L), marker, function(donor, at) {
    (donor * 7L + at * 3L + at %/% 2L + 1L) %% 2L
  })
  rownames(h1) <- rownames(h2) <- paste0("donor", seq_len(4L))
  colnames(h1) <- colnames(h2) <- paste0(chromosome, "_v", marker)
  args <- list(
    reference_haplotypes_h1 = h1,
    reference_haplotypes_h2 = h2,
    donor_population = c("A", "A", "B", "B"),
    ancestry_weights = c(A = 0.35, B = 0.65),
    N = c(A = 2, B = 2), Ne = c(A = 4, B = 7),
    rho = c(A = 0.8, B = 1.3),
    genetic_position = seq(0, 1.6, length.out = marker_count),
    mutation_age = rep(c(0.5, 2, 20, 1e6), length.out = marker_count),
    n = n, seed = seed, chromosome = rep(chromosome, marker_count),
    individual_offset = individual_offset
  )
  list(
    raw = do.call(gsim:::.gsim_hapnest_founders, c(args, list(
      return_genotypes = return_genotypes, return_segments = return_segments
    ))),
    packed = do.call(gsim:::.gsim_hapnest_founders_packed_chromosome,
                     c(list(backend = backend), args, list(
                       return_genotypes = return_genotypes,
                       return_segments = return_segments
                     )))
  )
}

.packed_pedigree_fixture <- function(backend, seed = 811, batch_size = NULL,
                                     chromosome = "chrM",
                                     return_genotypes = TRUE,
                                     return_crossovers = TRUE) {
  tab <- data.frame(
    animal = c("S", "D", "D2", "S2", "C", "FS", "PHS", "MHS", "G", "L"),
    sire = c(NA, NA, NA, NA, "S", "S", "S", "S2", "C", NA),
    dam = c(NA, NA, NA, NA, "D", "D", "D2", "D", "D2", NA),
    stringsAsFactors = FALSE
  )
  pedigree <- structure(list(
    pedigree = tab[c(10, 6, 3, 9, 1, 7, 4, 5, 2, 8), , drop = FALSE],
    canonical_order = tab$animal,
    external_order = tab$animal[c(10, 6, 3, 9, 1, 7, 4, 5, 2, 8)]
  ), class = "gsim_pedigree")
  marker_count <- 25L
  variant <- paste0(chromosome, "_v", seq_len(marker_count))
  h1 <- outer(seq_len(5L), seq_len(marker_count), function(i, j) {
    (i + j + j %/% 4L) %% 2L
  })
  h2 <- outer(seq_len(5L), seq_len(marker_count), function(i, j) {
    (i * 3L + j * 5L + j %/% 3L + 1L) %% 2L
  })
  rownames(h1) <- rownames(h2) <- c("S", "D", "D2", "S2", "L")
  colnames(h1) <- colnames(h2) <- variant
  founders <- list(h1 = h1, h2 = h2)
  position <- seq(0, 8, length.out = marker_count)
  args <- list(
    pedigree = pedigree, founder_haplotypes = founders,
    chromosome = rep(chromosome, marker_count),
    genetic_position = position, seed = seed,
    return_genotypes = return_genotypes,
    return_crossovers = return_crossovers, batch_size = batch_size
  )
  list(
    raw = do.call(gsim:::.gsim_pedigree_genotypes, args),
    packed = do.call(gsim:::.gsim_pedigree_genotypes_packed_chromosome,
                     c(list(backend = backend), args))
  )
}

.packed_mendelian_errors <- function(h1, h2, alignment) {
  errors <- 0L
  for (i in which(!alignment$founder)) {
    sire <- match(alignment$sire[[i]], alignment$animal)
    dam <- match(alignment$dam[[i]], alignment$animal)
    errors <- errors + sum(!(h1[i, ] == h1[sire, ] | h1[i, ] == h2[sire, ]))
    errors <- errors + sum(!(h2[i, ] == h1[dam, ] | h2[i, ] == h2[dam, ]))
  }
  errors
}

testthat::test_that("gbits backend fails clearly when unavailable", {
  testthat::expect_error(gsim:::.gsim_gbits_backend(""),
                         "shared-library path is required")
})

testthat::test_that("packed words follow marker-major LSB-first storage", {
  backend <- .packed_backend()
  values <- matrix(as.raw(0), 65L, 2L)
  values[c(1L, 64L, 65L), 1L] <- as.raw(1)
  packed <- gsim:::.gsim_gbits_pack(backend, values)
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(packed), values)
  testthat::expect_identical(
    as.integer(gsim:::.gsim_gbits_word(packed, 1L, 1L)),
    c(1L, 0L, 0L, 0L, 0L, 0L, 0L, 128L)
  )
  testthat::expect_identical(
    as.integer(gsim:::.gsim_gbits_word(packed, 1L, 2L)),
    c(1L, rep.int(0L, 7L))
  )
  testthat::expect_identical(unname(gsim:::.gsim_gbits_info(packed)),
                             c(65, 2, 2, 32))
})

testthat::test_that("packed founder output and event records exactly match raw", {
  backend <- .packed_backend()
  fixture <- .packed_founder_fixture(backend)
  packed_h1 <- gsim:::.gsim_gbits_unpack(fixture$packed$h1)
  packed_h2 <- gsim:::.gsim_gbits_unpack(fixture$packed$h2)
  testthat::expect_identical(packed_h1, fixture$raw$h1)
  testthat::expect_identical(packed_h2, fixture$raw$h2)
  testthat::expect_identical(fixture$packed$genotypes, fixture$raw$genotypes)
  testthat::expect_identical(fixture$packed$segments, fixture$raw$segments)
  testthat::expect_identical(fixture$packed$sample_ids, rownames(fixture$raw$h1))
  testthat::expect_identical(fixture$packed$variant_ids, colnames(fixture$raw$h1))
  testthat::expect_identical(
    as.integer(fixture$packed$genotypes),
    as.integer(packed_h1) + as.integer(packed_h2)
  )

  asymmetric_h1 <- matrix(0, 2L, 9L)
  asymmetric_h2 <- matrix(1, 2L, 9L)
  dimnames(asymmetric_h1) <- dimnames(asymmetric_h2) <- list(
    c("d1", "d2"), paste0("asym", seq_len(9L))
  )
  asymmetric <- gsim:::.gsim_hapnest_founders_packed_chromosome(
    backend, asymmetric_h1, asymmetric_h2, c("P", "P"), c(P = 1),
    c(P = 2), c(P = 2), c(P = 1), seq(0, 1, length.out = 9L),
    rep(1e12, 9L), 6L, 17, rep("phase", 9L)
  )
  testthat::expect_true(all(gsim:::.gsim_gbits_unpack(asymmetric$h1) == 0))
  testthat::expect_true(all(gsim:::.gsim_gbits_unpack(asymmetric$h2) == 1))
})

testthat::test_that("packed founder streams are reproducible and option invariant", {
  backend <- .packed_backend()
  first <- .packed_founder_fixture(backend, seed = 77)
  same <- .packed_founder_fixture(backend, seed = 77, return_genotypes = FALSE)
  no_audit <- .packed_founder_fixture(
    backend, seed = 77, return_genotypes = FALSE, return_segments = FALSE
  )
  different <- .packed_founder_fixture(backend, seed = 78)
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h1),
                             gsim:::.gsim_gbits_unpack(same$packed$h1))
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h2),
                             gsim:::.gsim_gbits_unpack(same$packed$h2))
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h1),
                             gsim:::.gsim_gbits_unpack(no_audit$packed$h1))
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h2),
                             gsim:::.gsim_gbits_unpack(no_audit$packed$h2))
  testthat::expect_null(same$packed$genotypes)
  testthat::expect_null(no_audit$packed$segments)
  testthat::expect_false(identical(
    gsim:::.gsim_gbits_unpack(first$packed$h1),
    gsim:::.gsim_gbits_unpack(different$packed$h1)
  ))

  batch1 <- .packed_founder_fixture(backend, n = 5L, seed = 77)
  batch2 <- .packed_founder_fixture(
    backend, n = 7L, seed = 77, individual_offset = 5L
  )
  testthat::expect_identical(
    rbind(gsim:::.gsim_gbits_unpack(batch1$packed$h1),
          gsim:::.gsim_gbits_unpack(batch2$packed$h1)),
    first$raw$h1
  )
  testthat::expect_identical(
    rbind(gsim:::.gsim_gbits_unpack(batch1$packed$h2),
          gsim:::.gsim_gbits_unpack(batch2$packed$h2)),
    first$raw$h2
  )
  testthat::expect_identical(
    rbind(batch1$packed$genotypes, batch2$packed$genotypes),
    first$raw$genotypes
  )
  segments <- rbind(batch1$packed$segments, batch2$packed$segments)
  rownames(segments) <- NULL
  expected <- first$raw$segments
  rownames(expected) <- NULL
  testthat::expect_identical(segments, expected)

  set.seed(103)
  expected_rng <- runif(4)
  set.seed(103)
  invisible(.packed_founder_fixture(backend, n = 2L, seed = 77))
  testthat::expect_identical(runif(4), expected_rng)
})

testthat::test_that("packed pedigree meiosis exactly matches raw across generations", {
  backend <- .packed_backend()
  fixture <- .packed_pedigree_fixture(backend)
  h1 <- gsim:::.gsim_gbits_unpack(fixture$packed$h1)
  h2 <- gsim:::.gsim_gbits_unpack(fixture$packed$h2)
  testthat::expect_identical(h1, fixture$raw$h1)
  testthat::expect_identical(h2, fixture$raw$h2)
  testthat::expect_identical(fixture$packed$genotypes, fixture$raw$genotypes)
  testthat::expect_identical(fixture$packed$crossover_audit,
                             fixture$raw$crossover_audit)
  testthat::expect_identical(fixture$packed$sample_ids, fixture$raw$sample_ids)
  testthat::expect_identical(fixture$packed$variant_ids, fixture$raw$variant_ids)
  testthat::expect_identical(fixture$packed$variant_map, fixture$raw$variant_map)
  testthat::expect_identical(fixture$packed$pedigree_alignment,
                             fixture$raw$pedigree_alignment)
  testthat::expect_identical(
    .packed_mendelian_errors(h1, h2, fixture$packed$pedigree_alignment), 0L
  )
})

testthat::test_that("fixed crossover materialization exactly matches the raw oracle", {
  backend <- .packed_backend()
  positions <- 0:5
  h1_values <- matrix(as.raw(c(0, 1, 0, 1, 0, 1)), 1L)
  h2_values <- matrix(as.raw(c(1, 0, 1, 0, 1, 0)), 1L)
  h1 <- gsim:::.gsim_gbits_pack(backend, h1_values)
  h2 <- gsim:::.gsim_gbits_pack(backend, h2_values)
  cases <- list(
    zero = numeric(), between = 2.5, exact = 3,
    multiple = c(1, 3, 4.5)
  )
  for (crossovers in cases) {
    destination <- gsim:::.gsim_gbits_zero(backend, 1L, 6L)
    boundaries <- vapply(crossovers, function(x) {
      which(positions >= x)[[1L]] - 1L
    }, integer(1L))
    gsim:::.gsim_gbits_make_gamete(
      destination, 1L, h1, h2, 1L, 1L, boundaries
    )
    expected <- gsim:::.gsim_meiosis_materialize(
      h1_values[1, ], h2_values[1, ], positions, crossovers, 1L
    )
    testthat::expect_identical(
      as.raw(gsim:::.gsim_gbits_unpack(destination)[1, ]), expected
    )
  }
})

testthat::test_that("packed pedigree streams are batch and option invariant", {
  backend <- .packed_backend()
  first <- .packed_pedigree_fixture(backend, seed = 909, batch_size = 1L)
  same <- .packed_pedigree_fixture(backend, seed = 909, batch_size = 99L,
                                   return_genotypes = FALSE)
  no_audit <- .packed_pedigree_fixture(
    backend, seed = 909, return_genotypes = FALSE, return_crossovers = FALSE
  )
  different <- .packed_pedigree_fixture(backend, seed = 910)
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h1),
                             gsim:::.gsim_gbits_unpack(same$packed$h1))
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h2),
                             gsim:::.gsim_gbits_unpack(same$packed$h2))
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h1),
                             gsim:::.gsim_gbits_unpack(no_audit$packed$h1))
  testthat::expect_identical(gsim:::.gsim_gbits_unpack(first$packed$h2),
                             gsim:::.gsim_gbits_unpack(no_audit$packed$h2))
  testthat::expect_identical(first$packed$crossover_audit,
                             same$packed$crossover_audit)
  testthat::expect_null(same$packed$genotypes)
  testthat::expect_null(no_audit$packed$crossover_audit)
  testthat::expect_false(identical(
    gsim:::.gsim_gbits_unpack(first$packed$h1),
    gsim:::.gsim_gbits_unpack(different$packed$h1)
  ))
  set.seed(81)
  expected_rng <- runif(5)
  set.seed(81)
  invisible(.packed_pedigree_fixture(backend, seed = 909))
  testthat::expect_identical(runif(5), expected_rng)
})

testthat::test_that("three chromosomes are exactly independent of orchestration order", {
  backend <- .packed_backend()
  labels <- c("chrZ", "01", "CHR1")
  forward <- lapply(labels, function(label) {
    .packed_founder_fixture(backend, n = 5L, seed = 606, chromosome = label)
  })
  reverse <- lapply(rev(labels), function(label) {
    .packed_founder_fixture(backend, n = 5L, seed = 606, chromosome = label)
  })
  names(forward) <- labels
  names(reverse) <- rev(labels)
  for (label in labels) {
    testthat::expect_identical(
      gsim:::.gsim_gbits_unpack(forward[[label]]$packed$h1),
      gsim:::.gsim_gbits_unpack(reverse[[label]]$packed$h1)
    )
    testthat::expect_identical(
      gsim:::.gsim_gbits_unpack(forward[[label]]$packed$h2),
      forward[[label]]$raw$h2
    )
    testthat::expect_identical(forward[[label]]$packed$segments,
                               reverse[[label]]$packed$segments)
  }

  pedigree_forward <- lapply(labels, function(label) {
    .packed_pedigree_fixture(backend, seed = 707, chromosome = label)
  })
  pedigree_reverse <- lapply(rev(labels), function(label) {
    .packed_pedigree_fixture(backend, seed = 707, chromosome = label)
  })
  names(pedigree_forward) <- labels
  names(pedigree_reverse) <- rev(labels)
  for (label in labels) {
    testthat::expect_identical(
      gsim:::.gsim_gbits_unpack(pedigree_forward[[label]]$packed$h1),
      pedigree_forward[[label]]$raw$h1
    )
    testthat::expect_identical(
      gsim:::.gsim_gbits_unpack(pedigree_forward[[label]]$packed$h2),
      gsim:::.gsim_gbits_unpack(pedigree_reverse[[label]]$packed$h2)
    )
    testthat::expect_identical(
      pedigree_forward[[label]]$packed$crossover_audit,
      pedigree_reverse[[label]]$packed$crossover_audit
    )
  }
})

testthat::test_that("bounded memory accounting shows one-bit payload reduction", {
  backend <- .packed_backend()
  individuals <- 128L
  markers <- 129L
  values <- matrix(as.raw(rep(c(0L, 1L), length.out = individuals * markers)),
                   individuals, markers)
  first <- gsim:::.gsim_gbits_pack(backend, values)
  second <- gsim:::.gsim_gbits_pack(backend, values)
  raw_bytes <- 2 * length(values)
  packed_bytes <- gsim:::.gsim_gbits_info(first)[[4L]] +
    gsim:::.gsim_gbits_info(second)[[4L]]
  testthat::expect_identical(raw_bytes, 33024)
  testthat::expect_identical(unname(packed_bytes), 4128)
  testthat::expect_identical(raw_bytes / packed_bytes, 8)

  ordinary <- .packed_founder_fixture(backend, return_genotypes = FALSE)
  testthat::expect_null(ordinary$packed$genotypes)
  testthat::expect_identical(
    ordinary$packed$memory$decoded_genotype_bytes, 0
  )
})
