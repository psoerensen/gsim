# Exact tests for the internal HAPNEST founder-haplotype core.

hapnest_fixture <- function(
  reference_h1 = matrix(c(
    0, 1, 0, 1,
    1, 0, 1, 0,
    0, 0, 1, 1,
    1, 1, 0, 0
  ), nrow = 4L, byrow = TRUE),
  reference_h2 = reference_h1,
  donor_population = c("A", "A", "B", "B"),
  ancestry_weights = c(A = 1, B = 1),
  N = c(A = 2, B = 2),
  Ne = c(A = 2, B = 2),
  rho = c(A = 1, B = 1),
  genetic_position = c(0, 0.1, 0.2, 0.3),
  mutation_age = rep(1e6, 4L),
  chromosome = rep("1", 4L),
  n = 3L,
  seed = 123,
  donor_phase = "hapnest",
  return_genotypes = TRUE,
  individual_offset = 0L
) {
  gsim:::.gsim_hapnest_founders(
    reference_haplotypes_h1 = reference_h1,
    reference_haplotypes_h2 = reference_h2,
    donor_population = donor_population,
    ancestry_weights = ancestry_weights,
    N = N,
    Ne = Ne,
    rho = rho,
    genetic_position = genetic_position,
    mutation_age = mutation_age,
    chromosome = chromosome,
    n = n,
    seed = seed,
    donor_phase = donor_phase,
    return_genotypes = return_genotypes,
    individual_offset = individual_offset
  )
}

testthat::test_that("distribution scales use Julia's shape-scale conventions", {
  scales <- gsim:::.gsim_hapnest_scales(N = 20, Ne = 100, T = 4, rho = 0.5)
  testthat::expect_identical(unname(scales), c(5, 0.25))
  testthat::expect_named(scales, c("gamma_scale", "exponential_scale"))
})

testthat::test_that("inclusive endpoint exactly follows HAPNEST overshoot rules", {
  endpoint <- gsim:::.gsim_hapnest_segment_endpoint
  map <- c(0, 0.1, 0.2, 0.3)

  testthat::expect_identical(endpoint(map, 1L, 4L, 0.01), 2L)
  testthat::expect_identical(endpoint(map, 1L, 4L, 0.1), 3L)
  testthat::expect_identical(endpoint(map, 1L, 4L, 100), 4L)
  testthat::expect_identical(endpoint(map, 4L, 4L, 0), 4L)
  testthat::expect_identical(endpoint(7, 1L, 1L, 1e-300), 1L)
  testthat::expect_identical(endpoint(map, 1L, 4L, .Machine$double.xmax), 4L)
})

testthat::test_that("fixed segment copying uses strict mutation-age filtering", {
  copied <- gsim:::.gsim_hapnest_copy_segment(
    donor = c(1, 1, 1, 0), mutation_age = c(11, 10, 9, 100),
    start = 1L, end = 4L, T = 10
  )
  testthat::expect_identical(as.integer(copied), c(1L, 0L, 0L, 0L))

  middle <- gsim:::.gsim_hapnest_copy_segment(
    donor = c(0, 1, 1, 0), mutation_age = rep(20, 4),
    start = 2L, end = 3L, T = 10
  )
  testthat::expect_identical(as.integer(middle), c(1L, 1L))
})

testthat::test_that("paired haplotypes produce exact homozygous and heterozygous counts", {
  h1 <- matrix(c(0, 1, 1, 0, 1, 0), nrow = 2L, byrow = TRUE)
  h2 <- matrix(c(0, 0, 1, 1, 1, 0), nrow = 2L, byrow = TRUE)
  genotype <- gsim:::.gsim_hapnest_pair(h1, h2)
  testthat::expect_type(genotype, "raw")
  testthat::expect_identical(
    matrix(as.integer(genotype), nrow = 2L),
    matrix(c(0L, 1L, 2L, 1L, 2L, 0L), nrow = 2L, byrow = TRUE)
  )
})

testthat::test_that("a one-variant chromosome is generated as one inclusive segment", {
  out <- hapnest_fixture(
    reference_h1 = matrix(1, nrow = 1L, ncol = 1L),
    reference_h2 = matrix(1, nrow = 1L, ncol = 1L),
    donor_population = "A", ancestry_weights = c(A = 1),
    N = c(A = 1), Ne = c(A = 1), rho = c(A = 1),
    genetic_position = 0, mutation_age = 1e9, chromosome = "7",
    n = 1L, seed = 8
  )
  testthat::expect_identical(as.integer(out$h1), 1L)
  testthat::expect_identical(as.integer(out$h2), 1L)
  testthat::expect_identical(as.integer(out$genotypes), 2L)
  testthat::expect_true(all(out$segments$start == 1L))
  testthat::expect_true(all(out$segments$end == 1L))
})

testthat::test_that("one-population generation retains phase and exact counts", {
  reference <- matrix(c(0, 1, 0, 1), nrow = 1L)
  out <- hapnest_fixture(
    reference_h1 = reference,
    reference_h2 = reference,
    donor_population = "A",
    ancestry_weights = c(A = 1),
    N = c(A = 1), Ne = c(A = 0.01), rho = c(A = 1),
    genetic_position = c(0, 0.1, 0.2, 0.3),
    mutation_age = rep(1e9, 4L), n = 2L
  )

  expected <- matrix(rep(c(0L, 1L, 0L, 1L), each = 2L), nrow = 2L)
  testthat::expect_type(out$h1, "raw")
  testthat::expect_identical(matrix(as.integer(out$h1), 2L), expected)
  testthat::expect_identical(matrix(as.integer(out$h2), 2L), expected)
  testthat::expect_identical(
    matrix(as.integer(out$genotypes), 2L), 2L * expected
  )
  testthat::expect_true(all(out$segments$donor_population == "A"))
})

testthat::test_that("asymmetric H1 and H2 panels use HAPNEST phase-specific lookup", {
  h1 <- matrix(c(0, 0, 0, 0), nrow = 1L)
  h2 <- matrix(c(1, 1, 1, 1), nrow = 1L)
  out <- hapnest_fixture(
    reference_h1 = h1,
    reference_h2 = h2,
    donor_population = "A",
    ancestry_weights = c(A = 1),
    N = c(A = 1), Ne = c(A = 1), rho = c(A = 1),
    genetic_position = c(0, 0.1, 0.2, 0.3),
    mutation_age = rep(.Machine$double.xmax, 4L),
    n = 3L,
    seed = 404
  )

  testthat::expect_true(all(as.integer(out$h1) == 0L))
  testthat::expect_true(all(as.integer(out$h2) == 1L))
  testthat::expect_true(all(as.integer(out$genotypes) == 1L))
  testthat::expect_identical(out$settings$donor_phase, "hapnest")
  testthat::expect_identical(
    out$settings$reference_layout,
    "paired H1/H2 rows are reference individuals"
  )
  testthat::expect_true(all(out$segments$donor_individual == 1L))
})

testthat::test_that("multiple donor populations obey ancestry weights", {
  reference <- rbind(rep(0, 5), rep(1, 5))
  out <- hapnest_fixture(
    reference_h1 = reference,
    reference_h2 = reference,
    donor_population = c("A", "B"),
    ancestry_weights = c(A = 0, B = 4),
    N = c(A = 1, B = 1), Ne = c(A = 1, B = 0.01),
    rho = c(A = 1, B = 1), genetic_position = seq(0, 0.4, 0.1),
    mutation_age = rep(1e9, 5L), chromosome = rep("1", 5L), n = 4L
  )
  testthat::expect_true(all(as.integer(out$h1) == 1L))
  testthat::expect_true(all(as.integer(out$h2) == 1L))
  testthat::expect_true(all(out$segments$donor_population == "B"))

  mixed <- hapnest_fixture(n = 20L, seed = 5)
  testthat::expect_setequal(unique(mixed$segments$donor_population), c("A", "B"))
})

testthat::test_that("chromosome blocks are independent and never crossed", {
  out <- hapnest_fixture(
    chromosome = c("1", "1", "2", "2"),
    genetic_position = c(0, 0.1, 0, 0.1),
    mutation_age = rep(1e9, 4L), n = 5L
  )
  chr <- c("1", "1", "2", "2")
  testthat::expect_true(all(chr[out$segments$start] == out$segments$chromosome))
  testthat::expect_true(all(chr[out$segments$end] == out$segments$chromosome))
  testthat::expect_true(all(out$segments$start <= out$segments$end))
})

testthat::test_that("fixed seeds are reproducible, batch-stable, and ignore R RNG", {
  first <- hapnest_fixture(n = 4L, seed = 999)
  second <- hapnest_fixture(n = 4L, seed = 999)
  testthat::expect_identical(first$h1, second$h1)
  testthat::expect_identical(first$h2, second$h2)
  testthat::expect_equal(first$segments, second$segments, tolerance = 0)

  batch1 <- hapnest_fixture(n = 2L, seed = 999, individual_offset = 0L)
  batch2 <- hapnest_fixture(n = 2L, seed = 999, individual_offset = 2L)
  testthat::expect_identical(rbind(batch1$h1, batch2$h1), first$h1)
  testthat::expect_identical(rbind(batch1$h2, batch2$h2), first$h2)
  combined_segments <- rbind(batch1$segments, batch2$segments)
  row.names(combined_segments) <- NULL
  expected_segments <- first$segments
  row.names(expected_segments) <- NULL
  testthat::expect_equal(combined_segments, expected_segments, tolerance = 0)

  set.seed(71)
  expected_random <- runif(3)
  set.seed(71)
  invisible(hapnest_fixture(n = 2L, seed = 999))
  testthat::expect_identical(runif(3), expected_random)
})

testthat::test_that("the SplitMix64 stream has a fixed known answer", {
  out <- hapnest_fixture(n = 2L, seed = 123)
  testthat::expect_identical(
    matrix(as.integer(out$h1), nrow = 2L),
    matrix(c(1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L), nrow = 2L)
  )
  testthat::expect_identical(
    matrix(as.integer(out$h2), nrow = 2L),
    matrix(c(0L, 1L, 0L, 1L, 1L, 0L, 1L, 1L), nrow = 2L)
  )
  testthat::expect_identical(out$segments$donor_individual,
                             c(4L, 3L, 4L, 4L, 1L))
  testthat::expect_identical(out$segments$start, c(1L, 1L, 1L, 1L, 3L))
  testthat::expect_identical(out$segments$end, c(4L, 4L, 4L, 2L, 4L))
  testthat::expect_equal(
    out$segments$coalescent_age,
    c(2.80906884665165624, 0.19192326337633564, 0.37243645910795636,
      2.29588737857102476, 4.89290165842939384),
    tolerance = 1e-15
  )
})

testthat::test_that("optional genotype counts are omitted without changing haplotypes", {
  with_counts <- hapnest_fixture(seed = 321, return_genotypes = TRUE)
  without_counts <- hapnest_fixture(seed = 321, return_genotypes = FALSE)
  testthat::expect_null(without_counts$genotypes)
  testthat::expect_identical(with_counts$h1, without_counts$h1)
  testthat::expect_identical(with_counts$h2, without_counts$h2)
  testthat::expect_identical(
    with_counts$genotypes,
    gsim:::.gsim_hapnest_pair(with_counts$h1, with_counts$h2)
  )
})

testthat::test_that("invalid founder inputs fail explicitly", {
  testthat::expect_error(
    hapnest_fixture(genetic_position = c(0, 0.2, 0.1, 0.3)),
    "nondecreasing"
  )
  testthat::expect_error(
    hapnest_fixture(genetic_position = c(0, NA, 0.2, 0.3)), "finite"
  )
  testthat::expect_error(
    hapnest_fixture(chromosome = c("1", "2", "1", "2")), "contiguous"
  )
  testthat::expect_error(
    hapnest_fixture(ancestry_weights = c(A = 1, B = -1)), "nonnegative"
  )
  testthat::expect_error(
    hapnest_fixture(ancestry_weights = c(A = 0, B = 0)), "positive sum"
  )
  testthat::expect_error(
    hapnest_fixture(ancestry_weights = c(1, 1)), "uniquely named"
  )
  bad_reference <- matrix(c(0, 1, 2, 0), nrow = 1L)
  testthat::expect_error(
    hapnest_fixture(
      reference_h1 = bad_reference, reference_h2 = matrix(0, 1L, 4L),
      donor_population = "A",
      ancestry_weights = c(A = 1), N = c(A = 1), Ne = c(A = 1),
      rho = c(A = 1)
    ),
    "0/1"
  )
  testthat::expect_error(
    hapnest_fixture(
      donor_population = rep("A", 4), ancestry_weights = c(A = 1, B = 1)
    ),
    "has no reference individual"
  )
  testthat::expect_error(hapnest_fixture(mutation_age = c(1, 2, NA, 4)),
                         "finite nonnegative")
  testthat::expect_error(hapnest_fixture(donor_phase = "pooled"),
                         "pooled sampling is not implemented")
  testthat::expect_error(hapnest_fixture(N = c(A = 3, B = 2)),
                         "number of reference individuals")
  testthat::expect_error(
    hapnest_fixture(reference_h2 = matrix(0, nrow = 3L, ncol = 4L)),
    "identical dimensions"
  )
  named_h1 <- matrix(0, nrow = 4L, ncol = 4L,
                     dimnames = list(paste0("donor", 1:4), paste0("v", 1:4)))
  named_h2 <- named_h1
  rownames(named_h2) <- rev(rownames(named_h2))
  testthat::expect_error(
    hapnest_fixture(reference_h1 = named_h1, reference_h2 = named_h2),
    "row names must be identical"
  )
  named_h2 <- named_h1
  colnames(named_h2) <- rev(colnames(named_h2))
  testthat::expect_error(
    hapnest_fixture(reference_h1 = named_h1, reference_h2 = named_h2),
    "column names must be identical"
  )
})
