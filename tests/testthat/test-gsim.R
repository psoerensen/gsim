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
  rng_default <- .Random.seed
  explicit_zero <- gsim(
    W = W, architecture = "bayesr", n_causal = 7L, nt = 2L,
    rg = 0.2, h2 = c(0.4, 0.6), seed = 90210,
    scale_effects = FALSE, return_genotypes = TRUE,
    a = 0, b = 0, c = 0,
    ld_score = "unused", annotation_score = "unused"
  )
  rng_zero <- .Random.seed
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
  for (field in c("component", "B", "B_causal", "G", "E", "Y",
                  "causal_rsids", "causal_probability",
                  "marker_multipliers")) {
    testthat::expect_identical(sim[[field]], explicit_zero[[field]])
  }
  testthat::expect_identical(rng_default, rng_zero)
  testthat::expect_identical(
    sim$settings$marker_multipliers$variance_model, "constant"
  )
})

testthat::test_that("derived variance weights follow the frozen formula", {
  set.seed(106)
  W <- matrix(rbinom(120 * 8, 2, 0.3), 120, 8)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  p <- stats::setNames(seq(0.1, 0.45, length.out = 8), colnames(W))
  r <- stats::setNames(seq(0.5, 4, length.out = 8), colnames(W))
  s <- stats::setNames(seq(0.8, 2.2, length.out = 8), colnames(W))
  common <- list(
    W = W, architecture = "bayesr", causal_probability = 1,
    seed = 31415, scale_effects = FALSE
  )
  unit <- do.call(gsim, common)
  derived <- do.call(gsim, c(common, list(
    maf = p[rev(names(p))], ld_score = r[rev(names(r))],
    annotation_score = s[rev(names(s))], a = -0.4, b = -1, c = 0.5
  )))
  expected <- exp(
    -0.4 * log(p * (1 - p)) - log(r) + 0.5 * log(s)
  )

  testthat::expect_equal(
    derived$marker_multipliers, expected, tolerance = 1e-15
  )
  testthat::expect_identical(unit$component, derived$component)
  testthat::expect_identical(
    unit$marker_probabilities, derived$marker_probabilities
  )
  testthat::expect_equal(
    derived$B,
    sweep(unit$B, 1L, sqrt(expected), "*"),
    tolerance = 1e-15
  )
  testthat::expect_identical(
    derived$settings$marker_multipliers$exponents,
    c(a = -0.4, b = -1, c = 0.5)
  )
  testthat::expect_identical(
    derived$settings$marker_multipliers$variance_model, "derived"
  )
  testthat::expect_identical(
    derived$settings$marker_multipliers$sources$maf$source,
    "explicit_input"
  )
  testthat::expect_identical(
    derived$settings$marker_multipliers$sources$ld_score$source,
    "explicit_input"
  )
  testthat::expect_identical(
    derived$settings$marker_multipliers$sources$annotation_score$source,
    "explicit_input"
  )
  testthat::expect_equal(
    derived$causal$variance_weight, unname(expected), tolerance = 1e-15
  )
  testthat::expect_true(all(derived$causal$causal_probability == 1))

  maf_only <- do.call(gsim, c(common, list(maf = p, a = 0.25)))
  ld_only <- do.call(gsim, c(common, list(ld_score = r, b = -0.5)))
  annotation_only <- do.call(
    gsim, c(common, list(annotation_score = s, c = 2))
  )
  testthat::expect_equal(
    maf_only$marker_multipliers, (p * (1 - p))^0.25,
    tolerance = 1e-15
  )
  testthat::expect_equal(
    ld_only$marker_multipliers, r^-0.5, tolerance = 1e-15
  )
  testthat::expect_equal(
    annotation_only$marker_multipliers, s^2, tolerance = 1e-15
  )
})

testthat::test_that("derived variance inputs align and reject invalid models", {
  W <- matrix(rep(0:2, length.out = 90 * 6), 90, 6)
  colnames(W) <- paste0("marker", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  ids <- colnames(W)
  p <- stats::setNames(seq(0.1, 0.4, length.out = 6), ids)
  r <- stats::setNames(seq(1, 3, length.out = 6), ids)
  s <- stats::setNames(seq(0.7, 1.7, length.out = 6), ids)
  common <- list(
    W = W, architecture = "bayesc", causal_probability = 1,
    seed = 2718, scale_effects = FALSE
  )
  canonical <- do.call(gsim, c(common, list(
    maf = p, ld_score = r, annotation_score = s, a = 0.2, b = -1, c = 1
  )))
  reordered <- do.call(gsim, c(common, list(
    maf = rev(p), ld_score = rev(r), annotation_score = rev(s),
    a = 0.2, b = -1, c = 1
  )))
  testthat::expect_identical(canonical$marker_multipliers,
                             reordered$marker_multipliers)
  testthat::expect_identical(canonical$B, reordered$B)

  testthat::expect_error(
    do.call(gsim, c(common, list(a = 1, marker_multipliers = p))),
    "cannot be combined"
  )
  testthat::expect_error(
    do.call(gsim, c(common, list(b = 1))), "ld_score is required"
  )
  testthat::expect_error(
    do.call(gsim, c(common, list(c = 1))), "annotation_score is required"
  )
  testthat::expect_error(
    do.call(gsim, c(common, list(maf = unname(p), a = 1))), "nonempty names"
  )
  bad_names <- p
  names(bad_names)[1L] <- "unexpected"
  testthat::expect_error(
    do.call(gsim, c(common, list(maf = bad_names, a = 1))),
    "missing: marker1.*unexpected: unexpected"
  )
  duplicated <- r
  names(duplicated)[1L] <- names(duplicated)[2L]
  testthat::expect_error(
    do.call(gsim, c(common, list(ld_score = duplicated, b = 1))),
    "duplicated: marker2"
  )
  for (bad in list(NA_real_, Inf, 0, -1)) {
    invalid <- s
    invalid[3L] <- bad
    testthat::expect_error(
      do.call(gsim, c(common, list(annotation_score = invalid, c = 1))),
      "marker3"
    )
  }
  invalid_maf <- p
  invalid_maf[4L] <- 1
  testthat::expect_error(
    do.call(gsim, c(common, list(maf = invalid_maf, a = 1))), "marker4"
  )
  testthat::expect_error(
    do.call(gsim, c(common, list(maf = p, a = 10000))),
    "invalid markers"
  )
  testthat::expect_error(
    do.call(gsim, c(common, list(a = c(1, 2)))), "finite numeric scalar"
  )

  scalar <- do.call(gsim, c(common, list(marker_multipliers = 2)))
  testthat::expect_equal(scalar$marker_multipliers,
                         stats::setNames(rep(2, 6), ids))
})

testthat::test_that("Glist metadata derives weights before causal-only loading", {
  set.seed(107)
  W_store <- matrix(rbinom(100 * 12, 2, 0.3), 100, 12)
  colnames(W_store) <- paste0("m", seq_len(ncol(W_store)))
  rownames(W_store) <- paste0("id", seq_len(nrow(W_store)))
  maf <- stats::setNames(seq(0.08, 0.41, length.out = 12), colnames(W_store))
  ld <- stats::setNames(seq(0.7, 3.45, length.out = 12), colnames(W_store))
  Glist <- list(
    ids = rownames(W_store),
    rsids = list(`1` = colnames(W_store)[1:6],
                 `2` = colnames(W_store)[7:12]),
    rsidsLD = list(`2` = colnames(W_store)[7:12],
                   `1` = colnames(W_store)[1:6]),
    maf = list(unname(maf[1:6]), unname(maf[7:12])),
    ldscores = list(ld[7:12], ld[1:6])
  )
  canonical_ids <- unname(unlist(Glist$rsidsLD, use.names = FALSE))
  q <- stats::setNames(rep(c(1, 0, 0), length.out = 12), canonical_ids)
  requested <- new.env(parent = emptyenv())
  requested$rsids <- character(0)
  fake_getG <- function(Glist, rsids, ids, chr = NULL,
                        impute = TRUE, scale = FALSE) {
    requested$rsids <- c(requested$rsids, rsids)
    W_store[ids, rsids, drop = FALSE]
  }
  common <- list(
    Glist = Glist, architecture = "bayesr", causal_probability = q,
    a = -0.4, b = -1, seed = 1618, scale_effects = FALSE,
    getG_fun = fake_getG
  )
  from_glist <- do.call(gsim, common)
  requested_from_glist <- requested$rsids
  requested$rsids <- character(0)
  explicit <- do.call(gsim, c(common, list(
    maf = maf[canonical_ids], ld_score = ld[canonical_ids]
  )))
  expected <- (maf[canonical_ids] * (1 - maf[canonical_ids]))^-0.4 *
    ld[canonical_ids]^-1

  testthat::expect_equal(from_glist$marker_multipliers, expected,
                         tolerance = 1e-15)
  testthat::expect_identical(from_glist$marker_multipliers,
                             explicit$marker_multipliers)
  testthat::expect_identical(from_glist$component, explicit$component)
  testthat::expect_identical(from_glist$B, explicit$B)
  testthat::expect_identical(from_glist$Y, explicit$Y)
  testthat::expect_setequal(requested_from_glist, names(q)[q == 1])
  testthat::expect_setequal(requested$rsids, names(q)[q == 1])
  testthat::expect_length(requested_from_glist, sum(q))
  testthat::expect_identical(
    from_glist$settings$marker_multipliers$sources$maf$source, "Glist"
  )
  testthat::expect_identical(
    from_glist$settings$marker_multipliers$sources$ld_score$source, "Glist"
  )

  missing_ld <- Glist
  missing_ld$ldscores <- NULL
  bad_common <- common
  bad_common$Glist <- missing_ld
  testthat::expect_error(
    do.call(gsim, bad_common),
    "Glist\\$ldscores"
  )
})

testthat::test_that("scalar annotation variance remains distinct from SBayesRC", {
  set.seed(108)
  W <- matrix(rbinom(130 * 10, 2, 0.27), 130, 10)
  colnames(W) <- paste0("m", seq_len(ncol(W)))
  rownames(W) <- paste0("id", seq_len(nrow(W)))
  A <- cbind(binary = rep(c(0, 1), each = 5),
             continuous = seq(-1, 1, length.out = 10))
  rownames(A) <- colnames(W)
  alpha <- matrix(c(0.7, 0, 0, 0.2, 0, 0), nrow = 2L, byrow = TRUE)
  score <- stats::setNames(seq(0.5, 1.4, length.out = 10), colnames(W))
  common <- list(
    W = W, A = A, architecture = "bayesr",
    annotation_model = "sbayesrc", alpha = alpha, n_causal = 5L,
    seed = 1414, scale_effects = FALSE
  )
  unit <- do.call(gsim, common)
  weighted <- do.call(gsim, c(common, list(annotation_score = score, c = 1)))

  testthat::expect_identical(unit$marker_probabilities,
                             weighted$marker_probabilities)
  testthat::expect_identical(unit$component, weighted$component)
  active <- unit$component != 1L
  testthat::expect_equal(
    weighted$B[active, , drop = FALSE],
    unit$B[active, , drop = FALSE] * sqrt(score[active]),
    tolerance = 1e-15
  )
  testthat::expect_error(
    do.call(gsim, c(common, list(c = 1))), "annotation_score is required"
  )
  testthat::expect_identical(weighted$A, unit$A)
  testthat::expect_identical(weighted$alpha, unit$alpha)
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
