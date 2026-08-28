#' Simulate a deterministic multigenerational pedigree workload
#'
#' `gsim_pedigree()` creates a linear-storage pedigree domain for pedigree-model
#' and sparse-solver validation. It never constructs a relationship matrix. The
#' returned pedigree table is in reproducibly arbitrary external order;
#' `canonical_order` is topological and `mapping` links both orders.
#'
#' A supplied `seed` calls [set.seed()] and the call then advances the caller's
#' RNG state. Without a seed, the current R RNG stream is used. `RNGkind()` and
#' the supplied seed are recorded.
#'
#' @param n_animals Optional total number of animals. When supplied, animals are
#'   distributed as evenly as possible over `n_generations`.
#' @param n_generations Number of generations.
#' @param animals_per_generation Generation size when `n_animals` is `NULL`.
#' @param founder_generations Number of initial all-founder generations.
#' @param sires_per_generation,dams_per_generation Maximum sizes of the
#'   sex-specific parent pools from each eligible generation. `NULL` selects a
#'   small scaling default.
#' @param overlapping_generation_probability Probability of choosing an
#'   available parent from two generations back rather than one.
#' @param unknown_sire_probability,unknown_dam_probability Parent-specific
#'   missingness probabilities.
#' @param new_founder_probability Probability that a later animal has neither
#'   recorded parent.
#' @param unphenotyped_probability Probability that an animal is marked as
#'   unavailable for observed-record generation.
#' @param permute_external_order Reproducibly permute the returned pedigree rows.
#' @param seed Optional nonnegative integer seed.
#'
#' @return An object of class `gsim_pedigree` containing the externally ordered
#'   pedigree, canonical and external ID vectors, an explicit order mapping,
#'   settings, diagnostics, and deterministic checksums.
#' @examples
#' ped <- gsim_pedigree(n_generations = 5, animals_per_generation = 40,
#'                      sires_per_generation = 6, dams_per_generation = 12,
#'                      seed = 101)
#' head(ped$pedigree)
#' @export
gsim_pedigree <- function(
  n_animals = NULL,
  n_generations = 5L,
  animals_per_generation = 200L,
  founder_generations = 2L,
  sires_per_generation = NULL,
  dams_per_generation = NULL,
  overlapping_generation_probability = 0.20,
  unknown_sire_probability = 0.02,
  unknown_dam_probability = 0.02,
  new_founder_probability = 0.01,
  unphenotyped_probability = 0.15,
  permute_external_order = TRUE,
  seed = NULL
) {
  n_generations <- .gsim_ped_integer(n_generations, "n_generations", 2L)
  animals_per_generation <- .gsim_ped_integer(
    animals_per_generation, "animals_per_generation", 2L)
  if (!is.null(n_animals)) {
    n_animals <- .gsim_ped_integer(n_animals, "n_animals", 2L)
    if (n_animals < 2L * n_generations) {
      .gsim_stop("n_animals must allow at least two animals per generation.")
    }
  }
  founder_generations <- .gsim_ped_integer(
    founder_generations, "founder_generations", 1L, n_generations - 1L)
  overlap <- .gsim_ped_probability(
    overlapping_generation_probability, "overlapping_generation_probability")
  unknown_sire_probability <- .gsim_ped_probability(
    unknown_sire_probability, "unknown_sire_probability")
  unknown_dam_probability <- .gsim_ped_probability(
    unknown_dam_probability, "unknown_dam_probability")
  new_founder_probability <- .gsim_ped_probability(
    new_founder_probability, "new_founder_probability")
  unphenotyped_probability <- .gsim_ped_probability(
    unphenotyped_probability, "unphenotyped_probability")
  if (!is.logical(permute_external_order) || length(permute_external_order) != 1L ||
      is.na(permute_external_order)) {
    .gsim_stop("permute_external_order must be TRUE or FALSE.")
  }
  seed <- .gsim_ped_seed(seed)
  if (!is.null(seed)) set.seed(seed)
  rng_kind <- RNGkind()

  sizes <- .gsim_ped_generation_sizes(n_animals, n_generations,
                                      animals_per_generation)
  n <- sum(sizes)
  if (is.null(sires_per_generation)) {
    sires_per_generation <- max(2L, as.integer(ceiling(sqrt(max(sizes)) / 2)))
  } else {
    sires_per_generation <- .gsim_ped_integer(
      sires_per_generation, "sires_per_generation", 1L)
  }
  if (is.null(dams_per_generation)) {
    dams_per_generation <- max(2L, as.integer(ceiling(sqrt(max(sizes)))))
  } else {
    dams_per_generation <- .gsim_ped_integer(
      dams_per_generation, "dams_per_generation", 1L)
  }

  generation <- rep.int(seq_len(n_generations), sizes)
  cohort <- generation
  width <- max(6L, nchar(as.character(n)))
  animal <- sprintf(paste0("A%0", width, "d"), seq_len(n))
  sex <- character(n)
  for (g in seq_len(n_generations)) {
    at <- which(generation == g)
    sex[at] <- sample(rep(c("M", "F"), length.out = length(at)),
                      length(at), replace = FALSE)
  }
  generation_members <- split(seq_len(n), generation)
  male_by_generation <- lapply(generation_members, function(at) at[sex[at] == "M"])
  female_by_generation <- lapply(generation_members, function(at) at[sex[at] == "F"])
  sire_index <- dam_index <- rep.int(NA_integer_, n)
  for (g in seq.int(founder_generations + 1L, n_generations)) {
    at <- which(generation == g)
    ns <- length(at)
    sire <- .gsim_ped_parent_draw(
      male_by_generation[[g - 1L]],
      if (g > 2L) male_by_generation[[g - 2L]] else integer(),
      g, ns, sires_per_generation, overlap)
    dam <- .gsim_ped_parent_draw(
      female_by_generation[[g - 1L]],
      if (g > 2L) female_by_generation[[g - 2L]] else integer(),
      g, ns, dams_per_generation, overlap)
    if (ns >= 2L) {
      sire[2L] <- sire[1L]
      dam[2L] <- dam[1L]
    }
    new_founder <- stats::runif(ns) < new_founder_probability
    if (new_founder_probability > 0 && ns >= 5L && !any(new_founder)) {
      new_founder[5L] <- TRUE
    }
    unknown_sire <- stats::runif(ns) < unknown_sire_probability
    unknown_dam <- stats::runif(ns) < unknown_dam_probability
    unknown_sire[new_founder] <- TRUE
    unknown_dam[new_founder] <- TRUE
    unknown_sire[seq_len(min(2L, ns))] <- FALSE
    unknown_dam[seq_len(min(2L, ns))] <- FALSE
    if (unknown_sire_probability > 0 && ns >= 3L && !any(unknown_sire & !new_founder)) {
      unknown_sire[3L] <- TRUE
    }
    if (unknown_dam_probability > 0 && ns >= 4L && !any(unknown_dam & !new_founder)) {
      unknown_dam[4L] <- TRUE
    }
    sire[unknown_sire] <- NA_integer_
    dam[unknown_dam] <- NA_integer_
    sire_index[at] <- sire
    dam_index[at] <- dam
  }
  phenotyped <- stats::runif(n) >= unphenotyped_probability
  if (unphenotyped_probability > 0 && all(phenotyped)) phenotyped[n] <- FALSE
  if (!any(phenotyped)) phenotyped[1L] <- TRUE
  canonical_table <- data.frame(
    animal = animal,
    sire = ifelse(is.na(sire_index), NA_character_, animal[sire_index]),
    dam = ifelse(is.na(dam_index), NA_character_, animal[dam_index]),
    sex = sex,
    generation = as.integer(generation),
    cohort = as.integer(cohort),
    phenotyped = phenotyped,
    stringsAsFactors = FALSE
  )
  external_index <- if (permute_external_order) sample.int(n) else seq_len(n)
  if (permute_external_order && identical(external_index, seq_len(n))) {
    external_index <- c(external_index[-1L], external_index[1L])
  }
  external_order <- animal[external_index]
  mapping <- data.frame(
    animal = animal,
    canonical_index = seq_len(n),
    external_index = match(animal, external_order),
    stringsAsFactors = FALSE
  )
  sire_counts <- table(canonical_table$sire, useNA = "no")
  dam_counts <- table(canonical_table$dam, useNA = "no")
  pair <- paste(canonical_table$sire, canonical_table$dam, sep = "|")
  pair <- pair[!is.na(canonical_table$sire) & !is.na(canonical_table$dam)]
  later <- canonical_table$generation > founder_generations
  diagnostics <- list(
    animals = n,
    founders = sum(is.na(canonical_table$sire) & is.na(canonical_table$dam)),
    later_new_founders = sum(later & is.na(canonical_table$sire) &
                               is.na(canonical_table$dam)),
    generations = n_generations,
    generation_sizes = stats::setNames(sizes, seq_len(n_generations)),
    sires_used = length(unique(canonical_table$sire[!is.na(canonical_table$sire)])),
    dams_used = length(unique(canonical_table$dam[!is.na(canonical_table$dam)])),
    unknown_sires = sum(later & is.na(canonical_table$sire)),
    unknown_dams = sum(later & is.na(canonical_table$dam)),
    unphenotyped = sum(!canonical_table$phenotyped),
    paternal_half_sib_sires = sum(sire_counts >= 2L),
    maternal_half_sib_dams = sum(dam_counts >= 2L),
    full_sib_families = sum(table(pair) >= 2L),
    overlapping_parent_links = sum(
      !is.na(sire_index) & generation - generation[sire_index] == 2L,
      !is.na(dam_index) & generation - generation[dam_index] == 2L
    ),
    topological = all(c(sire_index, dam_index)[!is.na(c(sire_index, dam_index))] <
                        rep(seq_len(n), 2L)[!is.na(c(sire_index, dam_index))]),
    storage_bytes = as.numeric(utils::object.size(canonical_table))
  )
  settings <- list(
    n_animals = n,
    n_generations = n_generations,
    generation_sizes = sizes,
    founder_generations = founder_generations,
    sires_per_generation = sires_per_generation,
    dams_per_generation = dams_per_generation,
    overlapping_generation_probability = overlap,
    unknown_sire_probability = unknown_sire_probability,
    unknown_dam_probability = unknown_dam_probability,
    new_founder_probability = new_founder_probability,
    unphenotyped_probability = unphenotyped_probability,
    permute_external_order = permute_external_order,
    seed = seed,
    rng_kind = rng_kind,
    rng_state_policy = "supplied seed resets, then call advances caller RNG state"
  )
  checksums <- list(
    identifiers_and_parents = .gsim_checksum(
      `row.names<-`(canonical_table[c("animal", "sire", "dam")], NULL)),
    canonical_order = .gsim_checksum(animal),
    external_order = .gsim_checksum(external_order)
  )
  structure(list(
    pedigree = canonical_table[external_index, , drop = FALSE],
    canonical_order = animal,
    external_order = external_order,
    mapping = mapping,
    settings = settings,
    diagnostics = diagnostics,
    checksums = checksums
  ), class = "gsim_pedigree")
}

#' Simulate observation records from a gsim pedigree
#'
#' `gsim_pedigree_records()` creates one scalar-long model view from a supplied
#' [gsim_pedigree()] domain. Missing traits and times are represented by absent
#' rows, never imputed phenotypes. The fixed design is returned as sorted,
#' one-based sparse triplets rather than a dense matrix.
#'
#' Latent animal values are deterministic recursive solver-workload values.
#' They are deliberately not claimed to be exact draws from a numerator
#' relationship covariance: `truth$exact_relationship_draw` is `FALSE`, and
#' covariance recovery is not a validation target. They are appropriate for
#' comparing solvers supplied with identical equations and right-hand sides.
#'
#' A supplied `seed` calls [set.seed()] and the call then advances the caller's
#' RNG state. Optional longitudinal prediction records contain animal, time,
#' basis and fitted truth, but no residual or phenotype.
#'
#' @param pedigree A `gsim_pedigree` object.
#' @param model One of `"single_trait"`, `"multitrait"`, or `"longitudinal"`.
#' @param fixed_effects Optional finite vector with one value per generated
#'   design column, unnamed in design order or named by all design columns.
#' @param genetic_covariance Genetic covariance (single/multitrait) or random-
#'   regression coefficient covariance (longitudinal).
#' @param residual_covariance Residual covariance for single/multitrait models.
#' @param residual_variances Longitudinal residual-group variances.
#' @param random_regression_rank Basis rank from one to three.
#' @param basis_function Optional longitudinal function accepting the time vector
#'   and returning a finite matrix with `random_regression_rank` columns. `NULL`
#'   uses the documented shifted-Legendre basis.
#' @param minimum_records,maximum_records Longitudinal observed records per
#'   phenotyped animal.
#' @param residual_groups Positive number of residual groups.
#' @param prediction_records Whether to return longitudinal prediction records.
#' @param permute_record_order Reproducibly permute scalar record order.
#' @param seed Optional nonnegative integer seed.
#'
#' @return An object of class `gsim_pedigree_records` with the original pedigree,
#'   scalar records, sparse fixed-design triplets, optional basis and prediction
#'   records, model truth, settings, diagnostics, and deterministic checksums.
#' @examples
#' ped <- gsim_pedigree(n_generations = 4, animals_per_generation = 20, seed = 1)
#' dat <- gsim_pedigree_records(ped, model = "longitudinal",
#'                              prediction_records = TRUE, seed = 2)
#' head(dat$records)
#' @export
gsim_pedigree_records <- function(
  pedigree,
  model = c("single_trait", "multitrait", "longitudinal"),
  fixed_effects = NULL,
  genetic_covariance = NULL,
  residual_covariance = NULL,
  residual_variances = NULL,
  random_regression_rank = 3L,
  basis_function = NULL,
  minimum_records = 1L,
  maximum_records = 5L,
  residual_groups = 3L,
  prediction_records = FALSE,
  permute_record_order = TRUE,
  seed = NULL
) {
  model <- match.arg(model)
  if (!inherits(pedigree, "gsim_pedigree") ||
      !is.data.frame(pedigree$pedigree) || is.null(pedigree$canonical_order)) {
    .gsim_stop("pedigree must be a gsim_pedigree object.")
  }
  minimum_records <- .gsim_ped_integer(minimum_records, "minimum_records", 1L)
  maximum_records <- .gsim_ped_integer(maximum_records, "maximum_records",
                                       minimum_records)
  residual_groups <- .gsim_ped_integer(residual_groups, "residual_groups", 1L)
  random_regression_rank <- .gsim_ped_integer(
    random_regression_rank, "random_regression_rank", 1L, 3L)
  if (!is.logical(prediction_records) || length(prediction_records) != 1L ||
      is.na(prediction_records) || !is.logical(permute_record_order) ||
      length(permute_record_order) != 1L || is.na(permute_record_order)) {
    .gsim_stop("prediction_records and permute_record_order must be TRUE or FALSE.")
  }
  if (prediction_records && model != "longitudinal") {
    .gsim_stop("prediction_records are available only for the longitudinal model.")
  }
  if (!is.null(basis_function) && model != "longitudinal") {
    .gsim_stop("basis_function is available only for the longitudinal model.")
  }
  seed <- .gsim_ped_seed(seed)
  if (!is.null(seed)) set.seed(seed)
  rng_kind <- RNGkind()
  canonical <- pedigree$canonical_order
  animals <- pedigree$pedigree[match(canonical, pedigree$pedigree$animal), , drop = FALSE]
  if (anyNA(animals$animal) || anyDuplicated(animals$animal)) {
    .gsim_stop("pedigree canonical mapping is invalid.")
  }
  observed_animal <- which(animals$phenotyped)
  if (!length(observed_animal)) .gsim_stop("pedigree has no phenotyped animals.")
  n <- length(canonical)
  common_truth <- list(
    exact_relationship_draw = FALSE,
    purpose = "solver_workload",
    covariance_recovery_target = FALSE,
    latent_semantics = paste(
      "recursive parent-average values with independent conditional innovations;",
      "not an inbreeding-aware exact A covariance draw"
    )
  )

  if (model == "single_trait") {
    G <- .gsim_validate_covariance(genetic_covariance, 1L,
                                  "genetic_covariance", matrix(1, 1L, 1L))
    R <- .gsim_validate_covariance(residual_covariance, 1L,
                                  "residual_covariance", matrix(1.5, 1L, 1L))
    latent <- .gsim_ped_latent(pedigree, G)
    number <- length(observed_animal)
    records <- data.frame(
      record = sprintf("R%08d", seq_len(number)),
      observation_unit = sprintf("OU%08d", seq_len(number)),
      animal = canonical[observed_animal],
      trait = "trait1",
      value = NA_real_,
      residual_group = paste0("RG", sample.int(residual_groups, number, TRUE)),
      record_order = seq_len(number),
      time = NA_real_,
      basis_row = NA_integer_,
      stringsAsFactors = FALSE
    )
    permutation <- .gsim_record_permutation(number, permute_record_order)
    records <- records[permutation, , drop = FALSE]
    records$record_order <- seq_len(number)
    design <- .gsim_fixed_design(records, pedigree, 1L, fixed_effects)
    group_scale <- seq(0.75, 1.25, length.out = residual_groups)
    group_index <- as.integer(sub("RG", "", records$residual_group, fixed = TRUE))
    residual <- stats::rnorm(number, sd = sqrt(R[1L, 1L] * group_scale[group_index]))
    genetic <- latent[match(records$animal, canonical), 1L]
    phenotype <- design$contribution + genetic + residual
    records$value <- phenotype
    truth <- c(common_truth, list(
      fixed_effects = design$beta,
      latent_animal_values = latent,
      genetic_record_values = genetic,
      fixed_contribution = design$contribution,
      total_fitted_truth = design$contribution + genetic,
      residual_values = residual,
      phenotypes = phenotype,
      genetic_covariance = G,
      residual_covariance = R,
      residual_group_scales = group_scale
    ))
    basis <- NULL
    prediction <- NULL
    model_diagnostics <- list(
      observation_units = number,
      scalar_records = number,
      observation_patterns = c(trait1 = number),
      records_per_animal = table(factor(records$animal, levels = canonical)),
      residual_group_counts = table(records$residual_group)
    )
  } else if (model == "multitrait") {
    G_default <- matrix(c(1, 0.4, 0.4, 0.7), 2L, 2L)
    R_default <- matrix(c(0.8, 0.2, 0.2, 0.5), 2L, 2L)
    G <- .gsim_validate_covariance(genetic_covariance, 2L,
                                  "genetic_covariance", G_default)
    R <- .gsim_validate_covariance(residual_covariance, 2L,
                                  "residual_covariance", R_default)
    latent <- .gsim_ped_latent(pedigree, G)
    colnames(latent) <- c("trait1", "trait2")
    number_units <- length(observed_animal)
    if (number_units < 3L) {
      .gsim_stop("multitrait simulation requires at least three phenotyped animals.")
    }
    pattern <- sample(c("both", "trait1_only", "trait2_only"), number_units,
                      replace = TRUE, prob = c(0.55, 0.225, 0.225))
    pattern[seq_len(3L)] <- c("both", "trait1_only", "trait2_only")
    unit_group <- sample.int(residual_groups, number_units, TRUE)
    unit_id <- sprintf("OU%08d", seq_len(number_units))
    counts <- ifelse(pattern == "both", 2L, 1L)
    unit_row <- rep.int(seq_len(number_units), counts)
    trait <- unlist(lapply(pattern, function(x) {
      if (x == "both") c("trait1", "trait2") else sub("_only", "", x)
    }), use.names = FALSE)
    total <- length(trait)
    records <- data.frame(
      record = sprintf("R%08d", seq_len(total)),
      observation_unit = unit_id[unit_row],
      animal = canonical[observed_animal[unit_row]],
      trait = trait,
      value = NA_real_,
      residual_group = paste0("RG", unit_group[unit_row]),
      record_order = seq_len(total),
      time = NA_real_,
      basis_row = NA_integer_,
      stringsAsFactors = FALSE
    )
    group_scale <- seq(0.75, 1.25, length.out = residual_groups)
    joint_residual <- matrix(stats::rnorm(number_units * 2L), number_units, 2L) %*%
      chol(R)
    joint_residual <- joint_residual * sqrt(group_scale[unit_group])
    residual <- joint_residual[cbind(unit_row, match(trait, c("trait1", "trait2")))]
    permutation <- .gsim_record_permutation(total, permute_record_order)
    records <- records[permutation, , drop = FALSE]
    residual <- residual[permutation]
    records$record_order <- seq_len(total)
    design <- .gsim_fixed_design(records, pedigree, 2L, fixed_effects)
    genetic <- latent[cbind(match(records$animal, canonical),
                            match(records$trait, c("trait1", "trait2")))]
    phenotype <- design$contribution + genetic + residual
    records$value <- phenotype
    truth <- c(common_truth, list(
      fixed_effects = design$beta,
      latent_animal_values = latent,
      genetic_record_values = genetic,
      fixed_contribution = design$contribution,
      total_fitted_truth = design$contribution + genetic,
      residual_values = residual,
      phenotypes = phenotype,
      genetic_covariance = G,
      residual_covariance = R,
      residual_group_scales = group_scale,
      observation_pattern = stats::setNames(pattern, unit_id),
      pattern_counts = table(factor(pattern,
                                    levels = c("both", "trait1_only", "trait2_only")))
    ))
    basis <- NULL
    prediction <- NULL
    model_diagnostics <- list(
      observation_units = number_units,
      scalar_records = total,
      observation_patterns = truth$pattern_counts,
      records_per_animal = table(factor(records$animal, levels = canonical)),
      residual_group_counts = table(records$residual_group)
    )
  } else {
    Q_default <- matrix(c(1, 0.35, -0.15,
                          0.35, 0.55, 0.12,
                          -0.15, 0.12, 0.30), 3L, 3L, byrow = TRUE)
    Q <- .gsim_validate_covariance(
      genetic_covariance, random_regression_rank, "genetic_covariance",
      Q_default[seq_len(random_regression_rank),
                seq_len(random_regression_rank), drop = FALSE])
    if (is.null(residual_variances)) {
      residual_variances <- if (residual_groups == 3L) c(0.4, 0.8, 1.2) else
        seq(0.4, 1.2, length.out = residual_groups)
    }
    if (!is.numeric(residual_variances) ||
        !(length(residual_variances) %in% c(1L, residual_groups)) ||
        any(!is.finite(residual_variances)) || any(residual_variances <= 0)) {
      .gsim_stop("residual_variances must contain one or residual_groups positive values.")
    }
    residual_variances <- rep(as.numeric(residual_variances),
                              length.out = residual_groups)
    latent <- .gsim_ped_latent(pedigree, Q)
    colnames(latent) <- paste0("coefficient", seq_len(random_regression_rank))
    number_units <- length(observed_animal)
    count <- sample(seq.int(minimum_records, maximum_records), number_units, TRUE)
    if (minimum_records < maximum_records && number_units >= 2L &&
        length(unique(count)) == 1L) {
      count[1:2] <- c(minimum_records, maximum_records)
    }
    total <- sum(count)
    unit_row <- rep.int(seq_len(number_units), count)
    time <- stats::runif(total)
    group <- sample.int(residual_groups, total, TRUE)
    records <- data.frame(
      record = sprintf("R%08d", seq_len(total)),
      observation_unit = sprintf("OU%08d", seq_len(total)),
      animal = canonical[observed_animal[unit_row]],
      trait = "trait1",
      value = NA_real_,
      residual_group = paste0("RG", group),
      record_order = seq_len(total),
      time = time,
      basis_row = seq_len(total),
      stringsAsFactors = FALSE
    )
    basis <- .gsim_basis(time, random_regression_rank, basis_function)
    residual <- stats::rnorm(total, sd = sqrt(residual_variances[group]))
    permutation <- .gsim_record_permutation(total, permute_record_order)
    records <- records[permutation, , drop = FALSE]
    basis <- basis[permutation, , drop = FALSE]
    residual <- residual[permutation]
    records$record_order <- seq_len(total)
    records$basis_row <- seq_len(total)
    design <- .gsim_fixed_design(records, pedigree, 1L, fixed_effects)
    coefficient_row <- latent[match(records$animal, canonical), , drop = FALSE]
    genetic <- rowSums(basis * coefficient_row)
    phenotype <- design$contribution + genetic + residual
    records$value <- phenotype
    prediction <- NULL
    prediction_basis <- NULL
    prediction_genetic <- NULL
    if (prediction_records) {
      prediction_time <- rep(c(0.125, 0.875), each = n)
      prediction_animal <- rep(canonical, times = 2L)
      prediction_basis <- .gsim_basis(prediction_time, random_regression_rank,
                                      basis_function)
      prediction_genetic <- rowSums(
        prediction_basis * latent[match(prediction_animal, canonical), , drop = FALSE])
      prediction_design <- .gsim_fixed_design(
        data.frame(animal = prediction_animal, trait = "trait1",
                   stringsAsFactors = FALSE),
        pedigree, 1L, design$beta)
      prediction <- data.frame(
        prediction_record = sprintf("P%08d", seq_len(2L * n)),
        animal = prediction_animal,
        time = prediction_time,
        basis_row = seq_len(2L * n),
        genetic_value = prediction_genetic,
        fixed_contribution = prediction_design$contribution,
        total_fitted_truth = prediction_genetic + prediction_design$contribution,
        stringsAsFactors = FALSE
      )
    }
    truth <- c(common_truth, list(
      fixed_effects = design$beta,
      latent_animal_coefficients = latent,
      genetic_record_values = genetic,
      fixed_contribution = design$contribution,
      total_fitted_truth = design$contribution + genetic,
      residual_values = residual,
      phenotypes = phenotype,
      coefficient_covariance = Q,
      residual_variances = residual_variances,
      prediction_basis = prediction_basis,
      prediction_genetic_values = prediction_genetic
    ))
    model_diagnostics <- list(
      observation_units = total,
      scalar_records = total,
      observation_patterns = c(longitudinal = total),
      records_per_animal = table(factor(records$animal, levels = canonical)),
      residual_group_counts = table(records$residual_group),
      time_range = range(records$time),
      record_count_distribution = table(count),
      basis_dimensions = dim(basis),
      basis_rank = qr(basis)$rank,
      mean_squared_basis_magnitude = mean(basis^2)
    )
  }
  checksums <- .gsim_record_checksums(pedigree, records, basis, design$design)
  diagnostics <- c(model_diagnostics, list(
    fixed_design = design$diagnostics,
    checksums = checksums,
    finite_truth = all(is.finite(records$value)) &&
      all(is.finite(truth$genetic_record_values)) &&
      all(is.finite(truth$residual_values)),
    exact_relationship_draw = FALSE
  ))
  out <- list(
    model = model,
    pedigree = pedigree,
    records = records,
    fixed_design = design$design,
    basis = basis,
    prediction_records = prediction,
    truth = truth,
    settings = list(
      seed = seed,
      rng_kind = rng_kind,
      rng_state_policy = "supplied seed resets, then call advances caller RNG state",
      residual_groups = residual_groups,
      permute_record_order = permute_record_order,
      minimum_records = minimum_records,
      maximum_records = maximum_records,
      random_regression_rank = random_regression_rank
    ),
    diagnostics = diagnostics
  )
  .gsim_finalize_records(out)
}
