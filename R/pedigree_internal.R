# Internal helpers for pedigree and record workload simulation.

.gsim_ped_integer <- function(x, name, lower = 0L, upper = .Machine$integer.max) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < lower || x > upper) {
    .gsim_stop(name, " must be one finite integer in [", lower, ", ", upper, "].")
  }
  as.integer(x)
}

.gsim_ped_probability <- function(x, name) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x < 0 || x > 1) {
    .gsim_stop(name, " must be one finite probability in [0, 1].")
  }
  as.numeric(x)
}

.gsim_ped_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  .gsim_ped_integer(seed, "seed", 0L, .Machine$integer.max)
}

.gsim_checksum <- function(x) {
  bytes <- as.integer(serialize(x, NULL, ascii = FALSE, version = 2L))
  modulus <- 2147483629
  value <- 104729
  if (length(bytes)) {
    for (byte in bytes) value <- (value * 131 + byte + 1) %% modulus
  }
  sprintf("%08x", as.integer(value))
}

.gsim_validate_covariance <- function(x, dimension, name, default) {
  if (is.null(x)) x <- default
  if (dimension == 1L && length(x) == 1L && is.numeric(x)) {
    x <- matrix(as.numeric(x), 1L, 1L)
  } else {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
  }
  if (!identical(dim(x), c(dimension, dimension)) || any(!is.finite(x)) ||
      max(abs(x - t(x))) > 1e-10) {
    .gsim_stop(name, " must be a finite symmetric ", dimension, " by ",
               dimension, " matrix.")
  }
  if (inherits(try(chol(x), silent = TRUE), "try-error")) {
    .gsim_stop(name, " must be positive definite.")
  }
  x
}

.gsim_ped_generation_sizes <- function(n_animals, n_generations,
                                       animals_per_generation) {
  if (is.null(n_animals)) {
    rep.int(animals_per_generation, n_generations)
  } else {
    base <- n_animals %/% n_generations
    remainder <- n_animals %% n_generations
    base + as.integer(seq_len(n_generations) <= remainder)
  }
}

.gsim_ped_parent_draw <- function(recent, older, generation, number, pool_size,
                                  overlap_probability) {
  restrict <- function(x) {
    if (!length(x)) return(x)
    sample(x, min(length(x), pool_size), replace = FALSE)
  }
  recent <- restrict(recent)
  older <- restrict(older)
  if (!length(recent) && !length(older)) {
    .gsim_stop("No eligible parents exist for generation ", generation, ".")
  }
  use_older <- length(older) > 0L &
    (length(recent) == 0L | stats::runif(number) < overlap_probability)
  out <- integer(number)
  if (any(!use_older)) out[!use_older] <- sample(recent, sum(!use_older), TRUE)
  if (any(use_older)) out[use_older] <- sample(older, sum(use_older), TRUE)
  if (overlap_probability > 0 && length(older) && !any(use_older)) {
    out[1L] <- sample(older, 1L)
  }
  out
}

.gsim_ped_latent <- function(pedigree, covariance) {
  canonical <- pedigree$canonical_order
  tab <- pedigree$pedigree[match(canonical, pedigree$pedigree$animal), , drop = FALSE]
  dimension <- nrow(covariance)
  innovations <- matrix(stats::rnorm(length(canonical) * dimension),
                        nrow = length(canonical), ncol = dimension) %*% chol(covariance)
  values <- innovations
  sire_position <- match(tab$sire, canonical)
  dam_position <- match(tab$dam, canonical)
  for (i in seq_along(canonical)) {
    sire <- sire_position[i]
    dam <- dam_position[i]
    known <- c(sire, dam)
    known <- known[!is.na(known)]
    if (length(known) == 1L) {
      values[i, ] <- 0.5 * values[known, ] + sqrt(0.75) * innovations[i, ]
    } else if (length(known) == 2L) {
      values[i, ] <- 0.5 * (values[known[1L], ] + values[known[2L], ]) +
        sqrt(0.5) * innovations[i, ]
    }
  }
  rownames(values) <- canonical
  values
}

.gsim_ped_management_group <- function(animal_index, n_animals) {
  number <- max(2L, min(20L, as.integer(ceiling(n_animals / 50))))
  paste0("H", 1L + (animal_index - 1L) %% number)
}

.gsim_fixed_design <- function(records, pedigree, n_traits = 1L,
                               fixed_effects = NULL) {
  map <- pedigree$pedigree[match(records$animal, pedigree$pedigree$animal), , drop = FALSE]
  animal_index <- match(records$animal, pedigree$canonical_order)
  herd <- .gsim_ped_management_group(animal_index, length(pedigree$canonical_order))
  cohorts <- sort(unique(pedigree$pedigree$cohort))
  herds <- sort(unique(.gsim_ped_management_group(
    seq_along(pedigree$canonical_order), length(pedigree$canonical_order))))
  trait_levels <- paste0("trait", seq_len(n_traits))
  block_names <- unlist(lapply(trait_levels, function(trait) {
    paste0(trait, ":", c("intercept", "sex", paste0("cohort=", cohorts),
                         paste0("herd=", herds)))
  }), use.names = FALSE)
  trait_index <- match(records$trait, trait_levels)
  if (anyNA(trait_index)) .gsim_stop("Record traits do not match the design traits.")
  term_count <- 2L + length(cohorts) + length(herds)
  offset <- (trait_index - 1L) * term_count
  rows <- rep(seq_len(nrow(records)), each = 4L)
  columns <- as.vector(rbind(
    offset + 1L,
    offset + 2L,
    offset + 2L + match(map$cohort, cohorts),
    offset + 2L + length(cohorts) + match(herd, herds)
  ))
  values <- as.vector(rbind(
    rep.int(1, nrow(records)),
    ifelse(map$sex == "F", 1, -1),
    rep.int(1, nrow(records)),
    rep.int(1, nrow(records))
  ))
  ordering <- order(rows, columns, method = "radix")
  design <- list(
    nrow = nrow(records),
    ncol = length(block_names),
    row = as.integer(rows[ordering]),
    column = as.integer(columns[ordering]),
    value = as.numeric(values[ordering]),
    column_names = block_names
  )
  if (anyDuplicated(paste(design$row, design$column, sep = ":"))) {
    .gsim_stop("Internal error: duplicate fixed-design triplets.")
  }
  if (is.null(fixed_effects)) {
    beta <- numeric(length(block_names))
    for (trait in seq_len(n_traits)) {
      at <- (trait - 1L) * term_count
      beta[at + 1L] <- 2 + 0.5 * (trait - 1L)
      beta[at + 2L] <- 0.2 * trait
      beta[at + 2L + seq_along(cohorts)] <-
        seq(-0.15, 0.15, length.out = length(cohorts))
      beta[at + 2L + length(cohorts) + seq_along(herds)] <-
        seq(-0.1, 0.1, length.out = length(herds))
    }
    names(beta) <- block_names
  } else {
    if (!is.numeric(fixed_effects) || any(!is.finite(fixed_effects))) {
      .gsim_stop("fixed_effects must be a finite numeric vector.")
    }
    if (!is.null(names(fixed_effects))) {
      location <- match(block_names, names(fixed_effects))
      if (anyNA(location) || anyDuplicated(names(fixed_effects))) {
        .gsim_stop("Named fixed_effects must uniquely contain every design column.")
      }
      beta <- as.numeric(fixed_effects[location])
    } else {
      if (length(fixed_effects) != length(block_names)) {
        .gsim_stop("fixed_effects must have one value per design column.")
      }
      beta <- as.numeric(fixed_effects)
    }
    names(beta) <- block_names
  }
  contribution <- numeric(nrow(records))
  for (i in seq_along(design$value)) {
    contribution[design$row[i]] <- contribution[design$row[i]] +
      design$value[i] * beta[design$column[i]]
  }
  list(design = design, beta = beta, contribution = contribution,
       herd = herd,
       diagnostics = list(
         dimensions = c(nrow = design$nrow, ncol = design$ncol),
         stored_entries = length(design$value),
         maximum_row_width = if (design$nrow) max(tabulate(design$row)) else 0L
       ))
}

.gsim_record_permutation <- function(number, permute) {
  if (!permute || number < 2L) return(seq_len(number))
  out <- sample.int(number)
  if (identical(out, seq_len(number))) out <- c(out[-1L], out[1L])
  out
}

.gsim_record_checksums <- function(pedigree, records, basis, design) {
  canonical_tab <- pedigree$pedigree[
    match(pedigree$canonical_order, pedigree$pedigree$animal), , drop = FALSE]
  rownames(canonical_tab) <- NULL
  list(
    pedigree = .gsim_checksum(canonical_tab[c("animal", "sire", "dam")]),
    observation_identifiers = .gsim_checksum(
      records[c("record", "observation_unit")]),
    animal_to_record = .gsim_checksum(records[c("record", "animal")]),
    trait_and_residual_group = .gsim_checksum(
      records[c("record", "trait", "residual_group")]),
    times = .gsim_checksum(records$time),
    basis = .gsim_checksum(basis),
    phenotypes = .gsim_checksum(records$value),
    fixed_design_triplets = .gsim_checksum(
      list(design$row, design$column, design$value, design$column_names))
  )
}

.gsim_basis <- function(time, rank, basis_function = NULL) {
  if (is.null(basis_function)) {
    z <- 2 * time - 1
    out <- cbind(1, z, 0.5 * (3 * z^2 - 1))[, seq_len(rank), drop = FALSE]
  } else {
    if (!is.function(basis_function)) {
      .gsim_stop("basis_function must be NULL or a function.")
    }
    out <- as.matrix(basis_function(time))
    storage.mode(out) <- "double"
    if (!identical(dim(out), c(length(time), rank)) || any(!is.finite(out))) {
      .gsim_stop("basis_function must return a finite length(time) by ", rank,
                 " numeric matrix.")
    }
  }
  colnames(out) <- paste0("basis", seq_len(rank))
  out
}

.gsim_finalize_records <- function(out) {
  out$diagnostics$object_sizes_bytes <- c(
    pedigree = as.numeric(utils::object.size(out$pedigree)),
    records = as.numeric(utils::object.size(out$records)),
    fixed_design = as.numeric(utils::object.size(out$fixed_design)),
    basis = as.numeric(utils::object.size(out$basis)),
    truth = as.numeric(utils::object.size(out$truth))
  )
  out$diagnostics$object_sizes_bytes <- c(
    out$diagnostics$object_sizes_bytes,
    total_returned = as.numeric(utils::object.size(out))
  )
  class(out) <- "gsim_pedigree_records"
  out
}
