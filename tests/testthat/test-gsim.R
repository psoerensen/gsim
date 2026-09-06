# Focused tests for the public gsim() contract.

.gsim_without_multiplier_additions <- function(x) {
  x$marker_multipliers <- NULL
  x$settings$marker_multipliers <- NULL
  x
}

testthat::test_that("in-memory simulation is reproducible and internally exact", {
  set.seed(11)
  W <- matrix(rbinom(600 * 80, 2, 0.3), 600, 80)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))

  x <- gsim(
    W = W,
    architecture = "bayesr",
    n_causal = 12L,
    h2 = 0.5,
    seed = 123,
    return_genotypes = TRUE
  )
  y <- gsim(
    W = W,
    architecture = "bayesr",
    n_causal = 12L,
    h2 = 0.5,
    seed = 123,
    return_genotypes = TRUE
  )

  testthat::expect_equal(x$B, y$B)
  testthat::expect_equal(x$Y, y$Y)
  testthat::expect_equal(nrow(x$causal), 12L)
  testthat::expect_equal(
    unname(x$G),
    unname(x$W_causal %*% x$B_causal),
    tolerance = 1e-12
  )
  testthat::expect_lt(x$exactness$max_y_minus_g_plus_e, 1e-12)
  testthat::expect_lt(abs(x$vg_observed - 1), 1e-10)
})

testthat::test_that("the default seeded component and effect draws are unchanged", {
  W <- matrix(
    (seq_len(120 * 24) * 17L) %% 3L, 120, 24,
    dimnames = list(paste0("id", 1:120), paste0("m", 1:24))
  )
  sim <- gsim(
    W = W, architecture = "bayesr", n_causal = 7L, nt = 2L,
    rg = 0.2, h2 = c(0.4, 0.6), seed = 90210,
    scale_effects = FALSE, return_genotypes = TRUE
  )
  expected <- matrix(c(
    -0.077974193640383160364, 0.038552830156511480597,
    -0.18304567539546776067, -0.054396810479884531719,
    -0.013787213342484894163, 0.28964250823060849749,
    0.056612880563344243623, 0.013425355401567792835,
    -0.034777784856836461980, -0.085832263786816700990,
    0.00015569042568410165132, -0.079326625829301181114,
    0.062282996563945917934, 0.21700905951862298204
  ), nrow = 7L)

  testthat::expect_identical(
    sim$causal_rsids, c("m2", "m3", "m6", "m9", "m10", "m23", "m24")
  )
  testthat::expect_equal(unname(sim$B_causal), expected, tolerance = 0)
})

testthat::test_that("mixed annotations produce valid SBayesRC probabilities", {
  set.seed(22)
  W <- matrix(rbinom(500 * 100, 2, 0.25), 500, 100)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  A <- cbind(
    binary = rep(c(0, 1), each = 50),
    continuous = seq(-2, 2, length.out = 100)
  )
  rownames(A) <- colnames(W)
  alpha <- matrix(
    c(1.2, 0, 0,
      0.4, 0, 0),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(colnames(A), paste0("stick", 1:3))
  )

  sim <- gsim(
    W = W,
    A = A,
    architecture = "bayesr",
    annotation_model = "sbayesrc",
    alpha = alpha,
    n_causal = 20L,
    seed = 456
  )

  P <- sim$marker_probabilities
  active <- 1 - P[, 1L]
  testthat::expect_equal(
    unname(rowSums(P)),
    rep(1, nrow(P)),
    tolerance = 1e-12
  )
  testthat::expect_true(all(P >= 0 & P <= 1))
  testthat::expect_gt(mean(active[A[, "binary"] == 1]),
                      mean(active[A[, "binary"] == 0]))
  testthat::expect_equal(sim$annotation_types,
                         c(binary = "binary", continuous = "continuous"))
})

testthat::test_that("unit marker multipliers are an exact RNG-neutral reduction", {
  set.seed(101)
  W <- matrix(rbinom(180 * 36, 2, 0.3), 180, 36)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  unit <- stats::setNames(rep(1, ncol(W)), colnames(W))
  common <- list(
    W = W, architecture = "bayesr", n_causal = 9L, h2 = 0.45,
    seed = 2468, return_genotypes = TRUE, compute_sumstats = TRUE
  )

  default <- do.call(gsim, common)
  rng_default <- .Random.seed
  supplied <- do.call(gsim, c(common, list(marker_multipliers = unit)))
  rng_supplied <- .Random.seed

  testthat::expect_identical(
    .gsim_without_multiplier_additions(default),
    .gsim_without_multiplier_additions(supplied)
  )
  testthat::expect_identical(rng_default, rng_supplied)
  testthat::expect_identical(default$marker_multipliers, unit)
  testthat::expect_identical(supplied$marker_multipliers, unit)
  testthat::expect_identical(default$settings$marker_multipliers$policy, "unit")
  testthat::expect_identical(
    supplied$settings$marker_multipliers$policy, "supplied"
  )
})

testthat::test_that("marker multipliers align once and reject invalid inputs", {
  set.seed(102)
  W <- matrix(rbinom(100 * 18, 2, 0.28), 100, 18)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  canonical <- stats::setNames(
    exp(seq(-0.5, 0.5, length.out = ncol(W))), colnames(W)
  )
  reordered <- canonical[rev(names(canonical))]
  common <- list(
    W = W, architecture = "bayesr", n_causal = 5L,
    seed = 975, scale_effects = FALSE
  )

  sim <- do.call(gsim, c(common, list(marker_multipliers = reordered)))
  testthat::expect_identical(sim$marker_multipliers, canonical)
  testthat::expect_identical(
    sim$settings$marker_multipliers$alignment, "canonical_marker_order"
  )

  duplicated <- unknown <- nonpositive <- nonfinite <- canonical
  names(duplicated)[1L] <- names(duplicated)[2L]
  names(unknown)[1L] <- "unknown_marker"
  nonpositive[1L] <- 0
  nonfinite[1L] <- Inf
  invalid <- list(
    unnamed = unname(canonical),
    missing = canonical[-1L],
    duplicated = duplicated,
    unknown = unknown,
    nonpositive = nonpositive,
    nonfinite = nonfinite
  )
  for (value in invalid) {
    testthat::expect_error(
      do.call(gsim, c(common, list(marker_multipliers = value))),
      "marker_multipliers"
    )
  }
})

testthat::test_that("marker multipliers change only active effect scales", {
  set.seed(103)
  W <- matrix(rbinom(160 * 30, 2, 0.31), 160, 30)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  q <- stats::setNames(exp(seq(-0.8, 0.8, length.out = ncol(W))), colnames(W))
  common <- list(
    W = W, architecture = "bayesr", nt = 2L, rg = 0.25,
    n_causal = 8L, seed = 8642, scale_effects = FALSE
  )

  unit <- do.call(gsim, common)
  rng_unit <- .Random.seed
  scaled <- do.call(gsim, c(common, list(marker_multipliers = q)))
  rng_scaled <- .Random.seed
  active <- unit$component != 1L
  expected <- sweep(
    unit$B[active, , drop = FALSE], 1L, sqrt(q[active]), "*"
  )

  testthat::expect_identical(unit$causal_rsids, scaled$causal_rsids)
  testthat::expect_identical(unit$component, scaled$component)
  testthat::expect_identical(
    unit$marker_probabilities, scaled$marker_probabilities
  )
  testthat::expect_equal(scaled$B[active, , drop = FALSE], expected,
                         tolerance = 1e-15)
  testthat::expect_true(all(scaled$B[!active, , drop = FALSE] == 0))
  testthat::expect_identical(rng_unit, rng_scaled)
  testthat::expect_identical(scaled$marker_multipliers, q)
  testthat::expect_equal(
    scaled$settings$marker_multipliers$geometric_mean,
    exp(mean(log(q))), tolerance = 1e-15
  )
  testthat::expect_identical(
    scaled$settings$marker_multipliers$n_markers, length(q)
  )
  testthat::expect_identical(scaled$settings$marker_multipliers$minimum, min(q))
  testthat::expect_identical(scaled$settings$marker_multipliers$maximum, max(q))
  testthat::expect_identical(scaled$settings$marker_multipliers$all_ones, FALSE)
})

testthat::test_that("membership and variance truth remain independent", {
  set.seed(104)
  W <- matrix(rbinom(140 * 24, 2, 0.27), 140, 24)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  A <- cbind(binary = rep(c(0, 1), each = 12),
             continuous = seq(-1, 1, length.out = 24))
  rownames(A) <- colnames(W)
  alpha <- matrix(
    c(0.8, 0, 0, 0.25, 0, 0), nrow = 2L, byrow = TRUE,
    dimnames = list(colnames(A), paste0("stick", 1:3))
  )
  q <- stats::setNames(exp(seq(-0.4, 0.4, length.out = 24)), colnames(W))
  common <- list(
    W = W, A = A, architecture = "bayesr",
    annotation_model = "sbayesrc", alpha = alpha,
    n_causal = 7L, seed = 1357, scale_effects = FALSE
  )

  inclusion_only <- do.call(gsim, common)
  combined <- do.call(gsim, c(common, list(marker_multipliers = q)))

  testthat::expect_identical(
    inclusion_only$marker_probabilities, combined$marker_probabilities
  )
  testthat::expect_identical(
    inclusion_only$continuation_probabilities,
    combined$continuation_probabilities
  )
  testthat::expect_identical(inclusion_only$component, combined$component)
  testthat::expect_identical(combined$marker_multipliers, q)
  testthat::expect_identical(
    combined$settings$annotation_model, "sbayesrc"
  )
  testthat::expect_identical(
    combined$settings$marker_multipliers$policy, "supplied"
  )
})

testthat::test_that("direct causal probabilities and variance weights are independent", {
  set.seed(105)
  W_store <- matrix(rbinom(90 * 12, 2, 0.3), 90, 12)
  colnames(W_store) <- paste0("m", seq_len(ncol(W_store)))
  rownames(W_store) <- paste0("id", seq_len(nrow(W_store)))
  Glist <- list(
    ids = rownames(W_store),
    rsids = list(`2` = colnames(W_store)[7:12],
                 `1` = colnames(W_store)[1:6])
  )
  q <- stats::setNames(
    c(0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0), colnames(W_store)
  )
  w <- stats::setNames(
    exp(seq(-0.6, 0.6, length.out = ncol(W_store))), colnames(W_store)
  )
  requested <- new.env(parent = emptyenv())
  requested$rsids <- character(0)
  fake_getG <- function(Glist, rsids, ids, chr = NULL,
                        impute = TRUE, scale = FALSE) {
    requested$rsids <- c(requested$rsids, rsids)
    W_store[ids, rsids, drop = FALSE]
  }
  common <- list(
    Glist = Glist, architecture = "bayesr", seed = 24601,
    causal_probability = q[rev(names(q))], getG_fun = fake_getG,
    scale_effects = FALSE, return_genotypes = TRUE
  )

  unit <- do.call(gsim, common)
  rng_unit <- .Random.seed
  requested_unit <- requested$rsids
  requested$rsids <- character(0)
  weighted <- do.call(gsim, c(common, list(
    marker_multipliers = w[rev(names(w))]
  )))
  rng_weighted <- .Random.seed

  canonical_ids <- unname(unlist(Glist$rsids, use.names = FALSE))
  q_canonical <- q[canonical_ids]
  w_canonical <- w[canonical_ids]
  active <- q_canonical == 1
  expected_probability <- cbind(
    1 - q_canonical,
    q_canonical * 0.8,
    q_canonical * 0.18,
    q_canonical * 0.02
  )
  dimnames(expected_probability) <- dimnames(weighted$marker_probabilities)
  testthat::expect_equal(
    weighted$marker_probabilities, expected_probability, tolerance = 1e-15
  )
  testthat::expect_identical(weighted$causal_probability, q_canonical)
  testthat::expect_identical(
    weighted$causal_rsids, names(q_canonical)[active]
  )
  testthat::expect_setequal(requested_unit, names(q_canonical)[active])
  testthat::expect_setequal(requested$rsids, names(q_canonical)[active])
  testthat::expect_length(requested_unit, sum(active))
  testthat::expect_length(requested$rsids, sum(active))
  testthat::expect_identical(unit$component, weighted$component)
  testthat::expect_equal(
    weighted$B[active, , drop = FALSE],
    sweep(unit$B[active, , drop = FALSE], 1L,
          sqrt(w_canonical[active]), "*"),
    tolerance = 1e-15
  )
  testthat::expect_true(all(weighted$B[!active, , drop = FALSE] == 0))
  testthat::expect_identical(rng_unit, rng_weighted)
  testthat::expect_identical(
    weighted$settings$causal_probability$policy, "supplied_bernoulli"
  )
  testthat::expect_equal(
    weighted$settings$causal_probability$expected_n_causal, sum(q)
  )
  testthat::expect_identical(
    weighted$causal$causal_probability, unname(q_canonical[active])
  )
  testthat::expect_identical(
    weighted$causal$variance_weight, unname(w_canonical[active])
  )

  all_causal <- gsim(
    W = W_store, architecture = "bayesc", causal_probability = 1,
    seed = 24602, scale_effects = FALSE
  )
  testthat::expect_true(all(all_causal$component == 2L))
  testthat::expect_true(all(all_causal$causal_probability == 1))
  testthat::expect_identical(
    all_causal$settings$causal_probability$policy, "supplied_bernoulli"
  )
})

testthat::test_that("direct causal probabilities reject ambiguous models", {
  W <- matrix(rep(0:2, length.out = 80 * 8), 80, 8)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  q <- stats::setNames(rep(0.5, ncol(W)), colnames(W))
  bad_duplicate <- bad_name <- q
  names(bad_duplicate)[1L] <- names(bad_duplicate)[2L]
  names(bad_name)[1L] <- "unknown"
  invalid <- list(
    unname(q), q[-1L], bad_duplicate, bad_name,
    replace(q, 1L, NA_real_), replace(q, 1L, -0.1),
    replace(q, 1L, 1.1), replace(q, seq_along(q), 0)
  )
  for (value in invalid) {
    testthat::expect_error(
      gsim(W = W, architecture = "bayesc", causal_probability = value,
           seed = 11),
      "causal_probability"
    )
  }

  testthat::expect_error(
    gsim(W = W, causal_probability = q, n_causal = 3L, seed = 11),
    "n_causal"
  )
  testthat::expect_error(
    gsim(W = W, architecture = "fixed", causal_probability = q,
         beta = rep(1, ncol(W)), seed = 11),
    "fixed"
  )
  testthat::expect_error(
    gsim(W = W, architecture = "clustered", causal_probability = q,
         block_id = rep(1:2, each = 4), seed = 11),
    "clustered"
  )
  testthat::expect_error(
    gsim(W = W, causal_probability = q,
         annotation_model = "sbayesrc", seed = 11),
    "sbayesrc"
  )
  testthat::expect_error(
    gsim(W = W, architecture = "bayesc", causal_probability = q,
         pi = c(1, 0), seed = 11),
    "positive mass"
  )
})

testthat::test_that("Glist phenotype generation requests causal markers only", {
  set.seed(33)
  W_store <- matrix(rbinom(120 * 60, 2, 0.3), 120, 60)
  colnames(W_store) <- paste0("m", seq_len(ncol(W_store)))
  rownames(W_store) <- paste0("id", seq_len(nrow(W_store)))
  Glist <- list(
    ids = rownames(W_store),
    rsids = list(`1` = colnames(W_store)[1:30],
                 `2` = colnames(W_store)[31:60])
  )
  q <- stats::setNames(
    seq(0.5, 1.5, length.out = ncol(W_store)), colnames(W_store)
  )

  requested <- new.env(parent = emptyenv())
  requested$rsids <- character(0)
  fake_getG <- function(Glist, rsids, ids, chr = NULL,
                        impute = TRUE, scale = FALSE) {
    requested$rsids <- c(requested$rsids, rsids)
    W_store[ids, rsids, drop = FALSE]
  }

  sim <- gsim(
    Glist = Glist,
    architecture = "bayesc",
    n_causal = 7L,
    h2 = 0.4,
    seed = 789,
    marker_multipliers = q[rev(names(q))],
    getG_fun = fake_getG,
    compute_sumstats = FALSE,
    return_genotypes = FALSE
  )

  testthat::expect_setequal(unique(requested$rsids), sim$causal_rsids)
  testthat::expect_equal(length(unique(requested$rsids)), 7L)
  testthat::expect_false(any(names(sim) == "W_causal"))
  testthat::expect_equal(dim(sim$Y), c(120L, 1L))
  testthat::expect_identical(sim$marker_multipliers, q)
})


testthat::test_that("gsim works with the qgg PLINK fixture", {
  testthat::skip_if_not_installed("qgg")

  bedfiles <- system.file(
    "extdata",
    paste0("sample_chr", 1:2, ".bed"),
    package = "qgg"
  )
  bimfiles <- sub("\\.bed$", ".bim", bedfiles)
  famfiles <- sub("\\.bed$", ".fam", bedfiles)

  Glist <- suppressMessages(
    qgg::gprep(
      study = "gsim-test",
      bedfiles = bedfiles,
      bimfiles = bimfiles,
      famfiles = famfiles
    )
  )

  sim <- gsim(
    Glist = Glist,
    architecture = "bayesr",
    n_causal = 10L,
    seed = 1,
    return_genotypes = TRUE
  )

  testthat::expect_equal(nrow(sim$Y), 489L)
  testthat::expect_equal(sim$settings$m, 2000L)
  testthat::expect_equal(sim$settings$n_causal, 10L)
  testthat::expect_equal(
    unname(sim$G),
    unname(sim$W_causal %*% sim$B_causal),
    tolerance = 1e-12
  )
})
