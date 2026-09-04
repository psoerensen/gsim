# Exact tests for the internal marker-level Mendelian inheritance oracle.

.marker_pedigree <- function(table, external_order = rev(table$animal)) {
  table$animal <- as.character(table$animal)
  table$sire <- as.character(table$sire)
  table$dam <- as.character(table$dam)
  if (!"sex" %in% names(table)) table$sex <- "U"
  if (!"generation" %in% names(table)) table$generation <- 1L
  if (!"cohort" %in% names(table)) table$cohort <- table$generation
  if (!"phenotyped" %in% names(table)) table$phenotyped <- TRUE
  canonical <- table$animal
  external_order <- as.character(external_order)
  structure(list(
    pedigree = table[match(external_order, canonical), , drop = FALSE],
    canonical_order = canonical,
    external_order = external_order,
    mapping = data.frame(
      animal = canonical,
      canonical_index = seq_along(canonical),
      external_index = match(canonical, external_order),
      stringsAsFactors = FALSE
    ),
    settings = list(), diagnostics = list(), checksums = list()
  ), class = "gsim_pedigree")
}

.three_generation_marker_pedigree <- function(external_order = NULL) {
  tab <- data.frame(
    animal = c("S", "D", "D2", "C", "G"),
    sire = c(NA, NA, NA, "S", "C"),
    dam = c(NA, NA, NA, "D", "D2"),
    generation = c(1L, 1L, 2L, 2L, 3L),
    stringsAsFactors = FALSE
  )
  if (is.null(external_order)) external_order <- c("G", "D", "C", "S", "D2")
  .marker_pedigree(tab, external_order)
}

.marker_founders <- function(variant_ids = paste0("v", 1:8)) {
  h1 <- rbind(
    D2 = rep(c(0, 1), length.out = length(variant_ids)),
    S = rep(c(0, 0, 1, 1), length.out = length(variant_ids)),
    D = rep(c(1, 0, 1, 0), length.out = length(variant_ids))
  )
  h2 <- rbind(
    D2 = rep(c(1, 0), length.out = length(variant_ids)),
    S = rep(c(1, 1, 0, 0), length.out = length(variant_ids)),
    D = rep(c(0, 1, 0, 1), length.out = length(variant_ids))
  )
  colnames(h1) <- colnames(h2) <- variant_ids
  list(h1 = h1, h2 = h2)
}

.marker_run <- function(
  pedigree = .three_generation_marker_pedigree(),
  founders = .marker_founders(),
  chromosome = rep(c("chr1", "chr2"), each = 4L),
  position = rep(c(0, 0.5, 1, 1.5), 2L),
  seed = 91,
  ...
) {
  gsim:::.gsim_pedigree_genotypes(
    pedigree = pedigree,
    founder_haplotypes = founders,
    chromosome = chromosome,
    genetic_position = position,
    seed = seed,
    ...
  )
}

testthat::test_that("fixed-event materialization follows the frozen boundary rule", {
  materialize <- gsim:::.gsim_meiosis_materialize
  h1 <- c(0, 0, 0, 0)
  h2 <- c(1, 1, 1, 1)
  position <- c(0, 1, 2, 3)

  testthat::expect_identical(
    as.integer(materialize(h1, h2, position, numeric(), 1L)),
    c(0L, 0L, 0L, 0L)
  )
  testthat::expect_identical(
    as.integer(materialize(h1, h2, position, numeric(), 2L)),
    c(1L, 1L, 1L, 1L)
  )
  testthat::expect_identical(
    as.integer(materialize(h1, h2, position, 1.5, 1L)),
    c(0L, 0L, 1L, 1L)
  )
  testthat::expect_identical(
    as.integer(materialize(h1, h2, position, 2, 1L)),
    c(0L, 0L, 1L, 1L)
  )
  testthat::expect_identical(
    as.integer(materialize(h1, h2, position, c(0.5, 2.5), 1L)),
    c(0L, 1L, 1L, 0L)
  )
  testthat::expect_identical(
    as.integer(materialize(h1, h2, position, 0, 1L)),
    c(1L, 1L, 1L, 1L)
  )

  leading_same_h1 <- c(0, 0, 0, 0)
  leading_same_h2 <- c(0, 1, 1, 1)
  testthat::expect_identical(
    as.integer(materialize(
      leading_same_h1, leading_same_h2, position, 0.5, 1L
    )),
    c(0L, 1L, 1L, 1L)
  )
  trailing_same_h2 <- c(1, 1, 1, 0)
  testthat::expect_identical(
    as.integer(materialize(h1, trailing_same_h2, position, 2.5, 1L)),
    c(0L, 0L, 0L, 0L)
  )
})

testthat::test_that("one-marker and zero-length chromosomes have zero crossovers", {
  ped <- .marker_pedigree(data.frame(
    animal = c("S", "D", "C"),
    sire = c(NA, NA, "S"), dam = c(NA, NA, "D"),
    stringsAsFactors = FALSE
  ))
  h1 <- matrix(c(0, 1), nrow = 2L,
               dimnames = list(c("D", "S"), "only"))
  h2 <- matrix(c(1, 0), nrow = 2L,
               dimnames = list(c("D", "S"), "only"))
  out <- gsim:::.gsim_pedigree_genotypes(
    ped, list(h1 = h1, h2 = h2), "chr7", 0.25, seed = 12
  )
  testthat::expect_true(all(out$crossover_audit$meioses$crossover_count == 0L))
  testthat::expect_equal(nrow(out$crossover_audit$crossovers), 0L)
  testthat::expect_true(as.integer(out$h1["C", ]) %in% c(0L, 1L))
  testthat::expect_true(as.integer(out$h2["C", ]) %in% c(0L, 1L))
  testthat::expect_identical(as.integer(out$genotypes),
                             as.integer(out$h1) + as.integer(out$h2))

  zero_span <- .marker_run(position = rep(c(0.4, 0.4, 0.4, 0.4), 2L))
  testthat::expect_true(
    all(zero_span$crossover_audit$meioses$crossover_count == 0L)
  )
})

testthat::test_that("founders align by ID and phases remain paternal and maternal", {
  out <- .marker_run()
  founders <- .marker_founders()
  founder_ids <- rownames(founders$h1)
  testthat::expect_identical(
    matrix(as.integer(out$h1[founder_ids, ]), nrow = length(founder_ids),
           dimnames = dimnames(founders$h1)),
    matrix(as.integer(founders$h1), nrow = length(founder_ids),
           dimnames = dimnames(founders$h1))
  )
  testthat::expect_identical(
    matrix(as.integer(out$h2[founder_ids, ]), nrow = length(founder_ids),
           dimnames = dimnames(founders$h2)),
    matrix(as.integer(founders$h2), nrow = length(founder_ids),
           dimnames = dimnames(founders$h2))
  )
  testthat::expect_identical(
    unname(out$settings$phase), c("paternal gamete", "maternal gamete")
  )
  testthat::expect_true(all(
    out$crossover_audit$meioses$parental_side[
      out$crossover_audit$meioses$child_phase == "H1"
    ] == "paternal"
  ))
  testthat::expect_true(all(
    out$crossover_audit$meioses$parental_side[
      out$crossover_audit$meioses$child_phase == "H2"
    ] == "maternal"
  ))
  testthat::expect_identical(out$sample_ids, out$pedigree_alignment$animal)
  testthat::expect_identical(out$variant_ids, out$variant_map$variant)
  testthat::expect_identical(out$variant_map$chromosome,
                             rep(c("chr1", "chr2"), each = 4L))
  testthat::expect_identical(out$variant_map$genetic_position,
                             rep(c(0, 0.5, 1, 1.5), 2L))
  founder_result_shape <- c(founders, list(
    genotypes = NULL, segments = data.frame(), settings = list()
  ))
  connected <- .marker_run(founders = founder_result_shape)
  testthat::expect_identical(out$h1, connected$h1)
  testthat::expect_identical(out$h2, connected$h2)
})

.mendelian_inconsistencies <- function(out) {
  alignment <- out$pedigree_alignment
  errors <- 0L
  for (i in which(!alignment$founder)) {
    sire <- match(alignment$sire[[i]], alignment$animal)
    dam <- match(alignment$dam[[i]], alignment$animal)
    paternal_ok <- out$h1[i, ] == out$h1[sire, ] |
      out$h1[i, ] == out$h2[sire, ]
    maternal_ok <- out$h2[i, ] == out$h1[dam, ] |
      out$h2[i, ] == out$h2[dam, ]
    errors <- errors + sum(!paternal_ok) + sum(!maternal_ok)
  }
  errors
}

testthat::test_that("children and grandchildren have zero Mendelian inconsistencies", {
  out <- .marker_run()
  testthat::expect_identical(.mendelian_inconsistencies(out), 0L)
  testthat::expect_true(all(as.integer(out$h1) %in% 0:1))
  testthat::expect_true(all(as.integer(out$h2) %in% 0:1))
  testthat::expect_true(all(as.integer(out$genotypes) %in% 0:2))
  testthat::expect_identical(
    as.integer(out$genotypes), as.integer(out$h1) + as.integer(out$h2)
  )
  testthat::expect_true(any(out$pedigree_alignment$animal == "G" &
                              !out$pedigree_alignment$founder))
})

testthat::test_that("stable streams are reproducible and invariant", {
  first <- .marker_run(seed = 811, batch_size = 1L)
  same <- .marker_run(seed = 811, batch_size = 99L)
  testthat::expect_identical(first$h1, same$h1)
  testthat::expect_identical(first$h2, same$h2)
  testthat::expect_identical(first$genotypes, same$genotypes)
  testthat::expect_identical(first$crossover_audit, same$crossover_audit)

  permuted <- .marker_run(
    pedigree = .three_generation_marker_pedigree(
      c("S", "C", "D2", "G", "D")
    ),
    seed = 811
  )
  testthat::expect_identical(first$h1, permuted$h1)
  testthat::expect_identical(first$h2, permuted$h2)
  testthat::expect_identical(first$genotypes, permuted$genotypes)
  testthat::expect_identical(first$crossover_audit, permuted$crossover_audit)

  no_genotypes <- .marker_run(seed = 811, return_genotypes = FALSE)
  testthat::expect_identical(first$h1, no_genotypes$h1)
  testthat::expect_identical(first$h2, no_genotypes$h2)
  testthat::expect_identical(first$crossover_audit,
                             no_genotypes$crossover_audit)
  testthat::expect_null(no_genotypes$genotypes)

  no_crossovers <- .marker_run(seed = 811, return_crossovers = FALSE)
  testthat::expect_identical(first$h1, no_crossovers$h1)
  testthat::expect_identical(first$h2, no_crossovers$h2)
  testthat::expect_identical(first$genotypes, no_crossovers$genotypes)
  testthat::expect_null(no_crossovers$crossover_audit)

  counts_only <- .marker_run(seed = 811, return_haplotypes = FALSE)
  testthat::expect_null(counts_only$h1)
  testthat::expect_null(counts_only$h2)
  testthat::expect_identical(first$genotypes, counts_only$genotypes)

  founders <- .marker_founders()
  one_chr <- gsim:::.gsim_pedigree_genotypes(
    .three_generation_marker_pedigree(),
    list(h1 = founders$h1[, 1:4], h2 = founders$h2[, 1:4]),
    chromosome = rep("chr1", 4L),
    genetic_position = c(0, 0.5, 1, 1.5), seed = 811
  )
  testthat::expect_identical(first$h1[, 1:4], one_chr$h1)
  testthat::expect_identical(first$h2[, 1:4], one_chr$h2)
  testthat::expect_identical(first$genotypes[, 1:4], one_chr$genotypes)

  extended_tab <- data.frame(
    animal = c("S", "D", "D2", "C", "G", "UNRELATED"),
    sire = c(NA, NA, NA, "S", "C", NA),
    dam = c(NA, NA, NA, "D", "D2", NA),
    generation = c(1L, 1L, 2L, 2L, 3L, 3L),
    stringsAsFactors = FALSE
  )
  extended_founders <- lapply(.marker_founders(), function(x) {
    rbind(x, UNRELATED = rep(0, ncol(x)))
  })
  extended <- .marker_run(
    pedigree = .marker_pedigree(extended_tab),
    founders = extended_founders,
    seed = 811
  )
  testthat::expect_identical(first$h1[first$sample_ids, ],
                             extended$h1[first$sample_ids, ])
  testthat::expect_identical(first$h2[first$sample_ids, ],
                             extended$h2[first$sample_ids, ])

  set.seed(17)
  expected <- runif(4)
  set.seed(17)
  invisible(.marker_run(seed = 811))
  testthat::expect_identical(runif(4), expected)
})

testthat::test_that("different seeds change inheritance on a long map", {
  variants <- paste0("long", 1:41)
  first <- .marker_run(
    founders = .marker_founders(variants),
    chromosome = rep("long_chr", 41L), position = seq(0, 20, length.out = 41L),
    seed = 101
  )
  second <- .marker_run(
    founders = .marker_founders(variants),
    chromosome = rep("long_chr", 41L), position = seq(0, 20, length.out = 41L),
    seed = 102
  )
  testthat::expect_false(identical(first$h1, second$h1))
  testthat::expect_false(identical(first$h2, second$h2))
  testthat::expect_false(identical(first$crossover_audit,
                                   second$crossover_audit))
})

testthat::test_that("crossovers remain within their chromosome blocks", {
  out <- .marker_run(seed = 41)
  x <- out$crossover_audit$crossovers
  block <- out$chromosome_blocks
  if (nrow(x)) {
    index <- match(x$chromosome, block$chromosome)
    lower <- out$variant_map$genetic_position[block$start_variant[index]]
    upper <- out$variant_map$genetic_position[block$end_variant[index]]
    testthat::expect_true(all(x$genetic_position > lower))
    testthat::expect_true(all(x$genetic_position < upper))
  }
  testthat::expect_true(all(
    out$crossover_audit$meioses$crossover_count ==
      ave(
        rep.int(1L, nrow(x)),
        interaction(x$animal, x$parental_side, x$chromosome, drop = TRUE),
        FUN = length
      )[match(
        interaction(
          out$crossover_audit$meioses$animal,
          out$crossover_audit$meioses$parental_side,
          out$crossover_audit$meioses$chromosome,
          drop = TRUE
        ),
        interaction(x$animal, x$parental_side, x$chromosome, drop = TRUE)
      )] | out$crossover_audit$meioses$crossover_count == 0L
  ))
})

testthat::test_that("invalid marker inheritance inputs fail explicitly", {
  ped <- .three_generation_marker_pedigree()
  founders <- .marker_founders()

  absent <- list(h1 = founders$h1[-1, ], h2 = founders$h2[-1, ])
  testthat::expect_error(.marker_run(founders = absent), "missing: D2")
  extra <- lapply(founders, function(x) rbind(x, X = x[1, ]))
  testthat::expect_error(.marker_run(founders = extra), "extra: X")
  duplicated <- founders
  rownames(duplicated$h1)[2] <- rownames(duplicated$h1)[1]
  rownames(duplicated$h2)[2] <- rownames(duplicated$h2)[1]
  testthat::expect_error(.marker_run(founders = duplicated), "unique animal IDs")
  nonbinary <- founders
  nonbinary$h1[1, 1] <- 2
  testthat::expect_error(.marker_run(founders = nonbinary), "0/1")
  wrong_dimension <- founders
  wrong_dimension$h2 <- wrong_dimension$h2[, -1, drop = FALSE]
  testthat::expect_error(.marker_run(founders = wrong_dimension),
                         "identical dimensions")
  wrong_rows <- founders
  rownames(wrong_rows$h2) <- rev(rownames(wrong_rows$h2))
  testthat::expect_error(.marker_run(founders = wrong_rows),
                         "row names must be identical")
  wrong_columns <- founders
  colnames(wrong_columns$h2) <- rev(colnames(wrong_columns$h2))
  testthat::expect_error(.marker_run(founders = wrong_columns),
                         "column names must be identical")
  duplicate_variant <- founders
  colnames(duplicate_variant$h1)[2] <- colnames(duplicate_variant$h1)[1]
  colnames(duplicate_variant$h2)[2] <- colnames(duplicate_variant$h2)[1]
  testthat::expect_error(.marker_run(founders = duplicate_variant),
                         "unique variant IDs")
  testthat::expect_error(.marker_run(position = c(0, 1, 0.5, 2, 0, 1, 2, 3)),
                         "nondecreasing")
  testthat::expect_error(.marker_run(position = c(0, 1, NA, 2, 0, 1, 2, 3)),
                         "finite")
  testthat::expect_error(.marker_run(chromosome = rep(c("a", "b", "a", "b"), 2L)),
                         "contiguous")

  partial_tab <- ped$pedigree
  partial_tab$dam[partial_tab$animal == "C"] <- NA
  partial_ped <- ped
  partial_ped$pedigree <- partial_tab
  testthat::expect_error(.marker_run(pedigree = partial_ped),
                         "Animal 'C'.*missing dam")

  bad_order <- ped
  bad_order$canonical_order <- c("S", "D", "D2", "G", "C")
  testthat::expect_error(.marker_run(pedigree = bad_order),
                         "parent-before-offspring.*G")

  testthat::expect_error(
    gsim:::.gsim_meiosis_materialize(c(0, 1), c(1, 0), c(0, 1), 1.1, 1L),
    "inside the chromosome interval"
  )
})

.relationship_pairs <- function(groups, tab) {
  children <- which(!is.na(tab$sire))
  po <- cbind(
    rep(children, each = 2L),
    as.vector(rbind(match(tab$sire[children], tab$animal),
                    match(tab$dam[children], tab$animal)))
  )
  within <- function(index) {
    if (length(index) < 2L) return(matrix(integer(), ncol = 2L))
    t(utils::combn(index, 2L))
  }
  fs <- do.call(rbind, lapply(groups, within))
  phs <- rbind(
    as.matrix(expand.grid(groups[[1L]], groups[[2L]])),
    as.matrix(expand.grid(groups[[3L]], groups[[4L]]))
  )
  mhs <- rbind(
    as.matrix(expand.grid(groups[[1L]], groups[[3L]])),
    as.matrix(expand.grid(groups[[2L]], groups[[4L]]))
  )
  founders <- which(is.na(tab$sire) & is.na(tab$dam))
  list(
    parent_offspring = po,
    full_siblings = fs,
    paternal_half_siblings = phs,
    maternal_half_siblings = mhs,
    unrelated_founders = t(utils::combn(founders, 2L))
  )
}

.relationship_experiment <- function() {
  family_size <- 8L
  founders <- c("S1", "S2", "D1", "D2")
  family_names <- list(
    paste0("F11_", seq_len(family_size)),
    paste0("F12_", seq_len(family_size)),
    paste0("F21_", seq_len(family_size)),
    paste0("F22_", seq_len(family_size))
  )
  children <- unlist(family_names, use.names = FALSE)
  tab <- data.frame(
    animal = c(founders, children),
    sire = c(rep(NA_character_, 4L),
             rep(c("S1", "S1", "S2", "S2"), each = family_size)),
    dam = c(rep(NA_character_, 4L),
            rep(c("D1", "D2", "D1", "D2"), each = family_size)),
    stringsAsFactors = FALSE
  )
  ped <- .marker_pedigree(tab)
  marker_count <- 256L
  variant <- paste0("independent_", seq_len(marker_count))
  set.seed(731)
  h1 <- matrix(sample.int(2L, 4L * marker_count, replace = TRUE) - 1L, 4L,
               dimnames = list(founders, variant))
  h2 <- matrix(sample.int(2L, 4L * marker_count, replace = TRUE) - 1L, 4L,
               dimnames = list(founders, variant))
  out <- gsim:::.gsim_pedigree_genotypes(
    ped, list(h1 = h1, h2 = h2), chromosome = paste0("chr", seq_len(marker_count)),
    genetic_position = rep(0, marker_count), seed = 20260904,
    return_crossovers = FALSE
  )
  genotype <- matrix(as.integer(out$genotypes), nrow = nrow(out$genotypes),
                     dimnames = dimnames(out$genotypes))
  pairs <- .relationship_pairs(
    lapply(family_names, function(x) match(x, tab$animal)), tab
  )
  pair_values <- function(markers, pair) {
    z <- genotype[, markers, drop = FALSE] - 1
    relationship <- 2 * tcrossprod(z) / length(markers)
    relationship[pair]
  }
  values <- lapply(pairs, function(pair) pair_values(seq_len(marker_count), pair))
  blocks <- split(seq_len(marker_count), rep(seq_len(16L), each = 16L))
  block_means <- lapply(pairs, function(pair) {
    vapply(blocks, function(markers) mean(pair_values(markers, pair)), numeric(1L))
  })
  means <- vapply(values, mean, numeric(1L))
  standard_errors <- vapply(block_means, stats::sd, numeric(1L)) / sqrt(length(blocks))
  variation <- vapply(values, stats::sd, numeric(1L))
  list(
    means = means,
    standard_errors = standard_errors,
    lower = means - 4 * standard_errors,
    upper = means + 4 * standard_errors,
    pair_sd = variation,
    expected = c(
      parent_offspring = 0.5,
      full_siblings = 0.5,
      paternal_half_siblings = 0.25,
      maternal_half_siblings = 0.25,
      unrelated_founders = 0
    )
  )
}

testthat::test_that("realized relationships agree within Monte Carlo uncertainty", {
  result <- .relationship_experiment()
  testthat::expect_true(all(result$expected >= result$lower &
                              result$expected <= result$upper))
  testthat::expect_gt(result$pair_sd[["full_siblings"]], 0)
  testthat::expect_gt(result$pair_sd[["paternal_half_siblings"]], 0)
  testthat::expect_gt(result$pair_sd[["maternal_half_siblings"]], 0)
})
