# Focused tests for deterministic pedigree and scalar-record workloads.

.ped_fixture <- function(seed = 101L) {
  gsim_pedigree(
    n_generations = 5L,
    animals_per_generation = 40L,
    founder_generations = 2L,
    sires_per_generation = 6L,
    dams_per_generation = 12L,
    overlapping_generation_probability = 0.25,
    unknown_sire_probability = 0.05,
    unknown_dam_probability = 0.05,
    new_founder_probability = 0.03,
    unphenotyped_probability = 0.15,
    permute_external_order = TRUE,
    seed = seed
  )
}

.expect_sparse_design <- function(x) {
  design <- x$fixed_design
  testthat::expect_equal(design$nrow, nrow(x$records))
  testthat::expect_length(design$row, length(design$column))
  testthat::expect_length(design$row, length(design$value))
  testthat::expect_true(all(is.finite(design$value)))
  testthat::expect_true(all(design$row >= 1L & design$row <= design$nrow))
  testthat::expect_true(all(design$column >= 1L & design$column <= design$ncol))
  testthat::expect_equal(anyDuplicated(paste(design$row, design$column)), 0L)
  testthat::expect_equal(order(design$row, design$column), seq_along(design$row))
  testthat::expect_true(all(tabulate(design$row, nbins = design$nrow) == 4L))
  intercept <- grepl(":intercept$", design$column_names[design$column])
  testthat::expect_equal(tabulate(design$row[intercept], nbins = design$nrow),
                         rep.int(1L, design$nrow))
}

.contains_forbidden_dense_object <- function(x) {
  found <- FALSE
  walk <- function(value, name = "") {
    if (is.matrix(value) && nrow(value) > 20L && nrow(value) == ncol(value)) {
      found <<- TRUE
    }
    if (is.list(value)) {
      child_names <- names(value)
      if (is.null(child_names)) child_names <- rep.int("", length(value))
      for (i in seq_along(value)) walk(value[[i]], child_names[i])
    }
  }
  walk(x)
  found
}

testthat::test_that("small pedigree is reproducible, valid, and relational", {
  x <- .ped_fixture()
  y <- .ped_fixture()
  z <- .ped_fixture(102L)
  testthat::expect_identical(x, y)
  testthat::expect_false(identical(x$checksums, z$checksums))
  testthat::expect_s3_class(x, "gsim_pedigree")
  testthat::expect_equal(nrow(x$pedigree), 200L)
  testthat::expect_false(anyNA(x$pedigree$animal))
  testthat::expect_false(any(!nzchar(x$pedigree$animal)))
  testthat::expect_equal(anyDuplicated(x$pedigree$animal), 0L)
  testthat::expect_setequal(x$canonical_order, x$external_order)
  testthat::expect_identical(x$pedigree$animal, x$external_order)
  testthat::expect_false(identical(x$canonical_order, x$external_order))
  testthat::expect_equal(x$mapping$external_index,
                         match(x$mapping$animal, x$external_order))
  parent <- c(x$pedigree$sire, x$pedigree$dam)
  testthat::expect_true(all(parent[!is.na(parent)] %in% x$canonical_order))
  canonical <- x$pedigree[match(x$canonical_order, x$pedigree$animal), ]
  child_position <- rep(seq_len(nrow(canonical)), 2L)
  parent_position <- match(c(canonical$sire, canonical$dam), x$canonical_order)
  known <- !is.na(parent_position)
  testthat::expect_true(all(parent_position[known] < child_position[known]))
  testthat::expect_true(x$diagnostics$topological)
  testthat::expect_setequal(canonical$sex, c("M", "F"))
  testthat::expect_gt(x$diagnostics$paternal_half_sib_sires, 0L)
  testthat::expect_gt(x$diagnostics$maternal_half_sib_dams, 0L)
  testthat::expect_gt(x$diagnostics$full_sib_families, 0L)
  testthat::expect_gt(x$diagnostics$overlapping_parent_links, 0L)
  testthat::expect_gt(x$diagnostics$unknown_sires, 0L)
  testthat::expect_gt(x$diagnostics$unknown_dams, 0L)
  testthat::expect_gt(x$diagnostics$later_new_founders, 0L)
  testthat::expect_gt(x$diagnostics$unphenotyped, 0L)
  testthat::expect_true(all(is.finite(unlist(x$diagnostics))))
  testthat::expect_lt(x$diagnostics$storage_bytes, 1000 * nrow(x$pedigree))
})

testthat::test_that("pedigree inputs reject invalid and impossible settings", {
  testthat::expect_error(gsim_pedigree(n_generations = 1L), "n_generations")
  testthat::expect_error(gsim_pedigree(n_animals = 9L, n_generations = 5L),
                         "at least two")
  testthat::expect_error(gsim_pedigree(founder_generations = 5L),
                         "founder_generations")
  testthat::expect_error(gsim_pedigree(unknown_sire_probability = 1.1),
                         "unknown_sire_probability")
  testthat::expect_error(gsim_pedigree(sires_per_generation = 0L),
                         "sires_per_generation")
  testthat::expect_error(gsim_pedigree(seed = -1L), "seed")
})

testthat::test_that("single-trait scalar records preserve truth and mappings", {
  ped <- .ped_fixture()
  x <- gsim_pedigree_records(ped, "single_trait", seed = 201L)
  y <- gsim_pedigree_records(ped, "single_trait", seed = 201L)
  testthat::expect_identical(x, y)
  testthat::expect_identical(x$pedigree, ped)
  testthat::expect_equal(anyDuplicated(x$records$record), 0L)
  testthat::expect_true(all(x$records$animal %in% ped$canonical_order))
  testthat::expect_equal(nrow(x$records), sum(ped$pedigree$phenotyped))
  testthat::expect_false(any(!ped$pedigree$phenotyped &
                               ped$pedigree$animal %in% x$records$animal))
  testthat::expect_equal(x$records$value,
                         unname(x$truth$fixed_contribution +
                                  x$truth$genetic_record_values +
                                  x$truth$residual_values))
  testthat::expect_true(all(is.finite(x$truth$latent_animal_values)))
  testthat::expect_equal(nrow(x$truth$latent_animal_values), nrow(ped$pedigree))
  testthat::expect_false(x$truth$exact_relationship_draw)
  testthat::expect_identical(x$truth$purpose, "solver_workload")
  .expect_sparse_design(x)
  testthat::expect_false(.contains_forbidden_dense_object(x))
})

testthat::test_that("multitrait records use absent rows for missing traits", {
  ped <- .ped_fixture()
  x <- gsim_pedigree_records(ped, "multitrait", seed = 202L)
  y <- gsim_pedigree_records(ped, "multitrait", seed = 202L)
  testthat::expect_identical(x, y)
  testthat::expect_setequal(names(x$truth$pattern_counts),
                            c("both", "trait1_only", "trait2_only"))
  testthat::expect_true(all(x$truth$pattern_counts > 0L))
  unit_width <- table(x$records$observation_unit)
  testthat::expect_setequal(as.integer(unit_width), c(1L, 2L))
  testthat::expect_false(anyNA(x$records$value))
  testthat::expect_equal(x$records$value,
                         x$truth$fixed_contribution +
                           x$truth$genetic_record_values +
                           x$truth$residual_values)
  testthat::expect_gt(min(eigen(x$truth$genetic_covariance,
                                symmetric = TRUE)$values), 0)
  testthat::expect_gt(min(eigen(x$truth$residual_covariance,
                                symmetric = TRUE)$values), 0)
  testthat::expect_true(all(is.finite(x$truth$latent_animal_values)))
  .expect_sparse_design(x)
  testthat::expect_false(.contains_forbidden_dense_object(x))
})

testthat::test_that("longitudinal basis, trajectories, and predictions align", {
  ped <- .ped_fixture()
  x <- gsim_pedigree_records(ped, "longitudinal",
                             prediction_records = TRUE, seed = 203L)
  y <- gsim_pedigree_records(ped, "longitudinal",
                             prediction_records = TRUE, seed = 203L)
  z <- gsim_pedigree_records(ped, "longitudinal",
                             prediction_records = TRUE, seed = 204L)
  testthat::expect_identical(x, y)
  testthat::expect_false(identical(x$diagnostics$checksums,
                                   z$diagnostics$checksums))
  testthat::expect_equal(x$records$basis_row, seq_len(nrow(x$records)))
  testthat::expect_equal(nrow(x$basis), nrow(x$records))
  testthat::expect_equal(x$basis[, 1L], rep.int(1, nrow(x$basis)))
  testthat::expect_equal(x$basis[, 2L], 2 * x$records$time - 1)
  testthat::expect_equal(x$diagnostics$basis_rank, 3L)
  coefficients <- x$truth$latent_animal_coefficients[
    match(x$records$animal, ped$canonical_order), , drop = FALSE]
  testthat::expect_equal(x$truth$genetic_record_values,
                         rowSums(x$basis * coefficients))
  testthat::expect_equal(x$records$value,
                         x$truth$total_fitted_truth + x$truth$residual_values)
  testthat::expect_true(all(table(x$records$animal) >= 1L &
                              table(x$records$animal) <= 5L))
  testthat::expect_gt(length(unique(as.integer(table(x$records$animal)))), 1L)
  prediction <- x$prediction_records
  prediction_coefficients <- x$truth$latent_animal_coefficients[
    match(prediction$animal, ped$canonical_order), , drop = FALSE]
  testthat::expect_equal(prediction$genetic_value,
                         rowSums(x$truth$prediction_basis * prediction_coefficients))
  testthat::expect_equal(prediction$total_fitted_truth,
                         prediction$genetic_value + prediction$fixed_contribution)
  testthat::expect_false(any(c("value", "residual") %in% names(prediction)))
  .expect_sparse_design(x)
  testthat::expect_false(.contains_forbidden_dense_object(x))
})

testthat::test_that("a supplied longitudinal basis is validated and aligned", {
  ped <- .ped_fixture()
  custom <- function(time) cbind(1, sin(pi * time))
  x <- gsim_pedigree_records(
    ped, "longitudinal", random_regression_rank = 2L,
    basis_function = custom, minimum_records = 1L, maximum_records = 2L,
    prediction_records = TRUE, seed = 205L
  )
  testthat::expect_equal(x$basis[, 2L], sin(pi * x$records$time))
  testthat::expect_equal(x$truth$prediction_basis[, 2L],
                         sin(pi * x$prediction_records$time))
  testthat::expect_error(
    gsim_pedigree_records(ped, "longitudinal", random_regression_rank = 2L,
                          basis_function = function(time) cbind(1, time, time^2),
                          seed = 1L),
    "length\\(time\\) by 2"
  )
  testthat::expect_error(
    gsim_pedigree_records(ped, "single_trait", basis_function = custom, seed = 1L),
    "only for the longitudinal"
  )
})
