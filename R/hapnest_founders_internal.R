# Internal founder-haplotype interface.  The model and provenance are frozen in
# docs/design/hapnest_founder_model.md; this is intentionally not exported.

#' @useDynLib gsim, .registration = TRUE
NULL

.gsim_hapnest_named_parameter <- function(x, populations, name) {
  if (!is.numeric(x) || is.null(names(x)) || anyNA(names(x)) ||
      any(!nzchar(names(x))) || anyDuplicated(names(x))) {
    .gsim_stop(name, " must be a uniquely named numeric vector.")
  }
  idx <- match(populations, names(x))
  if (anyNA(idx)) {
    .gsim_stop(name, " is missing a configured donor population.")
  }
  value <- as.double(x[idx])
  if (any(!is.finite(value)) || any(value <= 0)) {
    .gsim_stop(name, " must be finite and strictly positive.")
  }
  value
}

.gsim_hapnest_integer_scalar <- function(x, name, lower = 0) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x != floor(x) || x > .Machine$integer.max) {
    .gsim_stop(name, " must be an integer-valued scalar in [", lower,
               ", ", .Machine$integer.max, "].")
  }
  as.integer(x)
}

.gsim_hapnest_raw_matrix <- function(x, name) {
  if (!is.matrix(x) || !length(x) || any(dim(x) == 0L)) {
    .gsim_stop(name, " must be a nonempty matrix.")
  }
  if (!(is.raw(x) || is.logical(x) || is.integer(x) || is.double(x))) {
    .gsim_stop(name, " must contain 0/1 alleles.")
  }
  if (!is.raw(x) && (anyNA(x) || any(!is.finite(x)) || any(!(x %in% c(0, 1))))) {
    .gsim_stop(name, " must contain only nonmissing 0/1 alleles.")
  }
  if (is.raw(x) && any(!(as.integer(x) %in% c(0L, 1L)))) {
    .gsim_stop(name, " must contain only 0/1 alleles.")
  }
  matrix(as.raw(x), nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
}

.gsim_hapnest_founders <- function(
  reference_haplotypes_h1,
  reference_haplotypes_h2,
  donor_population,
  ancestry_weights,
  N,
  Ne,
  rho,
  genetic_position,
  mutation_age,
  n,
  seed,
  chromosome = rep.int("1", ncol(reference_haplotypes_h1)),
  donor_phase = "hapnest",
  return_genotypes = TRUE,
  individual_offset = 0L
) {
  reference_haplotypes_h1 <- .gsim_hapnest_raw_matrix(
    reference_haplotypes_h1, "reference_haplotypes_h1"
  )
  reference_haplotypes_h2 <- .gsim_hapnest_raw_matrix(
    reference_haplotypes_h2, "reference_haplotypes_h2"
  )
  if (!identical(dim(reference_haplotypes_h1), dim(reference_haplotypes_h2))) {
    .gsim_stop("reference_haplotypes_h1 and reference_haplotypes_h2 must have identical dimensions.")
  }
  if (!is.null(rownames(reference_haplotypes_h1)) ||
      !is.null(rownames(reference_haplotypes_h2))) {
    if (!identical(rownames(reference_haplotypes_h1),
                   rownames(reference_haplotypes_h2))) {
      .gsim_stop("H1 and H2 donor-individual row names must be identical when supplied.")
    }
  }
  if (!is.null(colnames(reference_haplotypes_h1)) ||
      !is.null(colnames(reference_haplotypes_h2))) {
    if (!identical(colnames(reference_haplotypes_h1),
                   colnames(reference_haplotypes_h2))) {
      .gsim_stop("H1 and H2 variant column names must be identical when supplied.")
    }
  }
  if (!is.character(donor_phase) || length(donor_phase) != 1L ||
      is.na(donor_phase) || donor_phase != "hapnest") {
    .gsim_stop("donor_phase must be 'hapnest'; pooled sampling is not implemented in this milestone.")
  }
  donor_population <- as.character(donor_population)
  if (length(donor_population) != nrow(reference_haplotypes_h1) ||
      anyNA(donor_population) || any(!nzchar(donor_population))) {
    .gsim_stop(
      "donor_population must provide one nonmissing label per reference individual."
    )
  }

  if (!is.numeric(ancestry_weights) || is.null(names(ancestry_weights)) ||
      anyNA(names(ancestry_weights)) || any(!nzchar(names(ancestry_weights))) ||
      anyDuplicated(names(ancestry_weights)) ||
      any(!is.finite(ancestry_weights)) || any(ancestry_weights < 0)) {
    .gsim_stop(
      "ancestry_weights must be uniquely named, finite, nonnegative, and have a positive sum."
    )
  }
  weight_sum <- sum(ancestry_weights)
  if (!is.finite(weight_sum) || weight_sum <= 0) {
    .gsim_stop(
      "ancestry_weights must be uniquely named, finite, nonnegative, and have a positive sum."
    )
  }
  active <- names(ancestry_weights)[ancestry_weights > 0]
  weights <- as.double(ancestry_weights[active])
  weights <- weights / sum(weights)
  donor_codes <- match(donor_population, active)
  for (i in seq_along(active)) {
    if (!any(donor_codes == i, na.rm = TRUE)) {
      .gsim_stop("Positive-weight donor population '", active[[i]],
                 "' has no reference individual.")
    }
  }
  donor_codes[is.na(donor_codes)] <- 0L

  N <- .gsim_hapnest_named_parameter(N, active, "N")
  Ne <- .gsim_hapnest_named_parameter(Ne, active, "Ne")
  rho <- .gsim_hapnest_named_parameter(rho, active, "rho")
  observed_N <- vapply(seq_along(active), function(i) {
    sum(donor_codes == i)
  }, integer(1L))
  if (!identical(N, as.double(observed_N))) {
    .gsim_stop("N must equal the number of reference individuals in each active donor population.")
  }

  m <- ncol(reference_haplotypes_h1)
  genetic_position <- as.double(genetic_position)
  mutation_age <- as.double(mutation_age)
  if (length(genetic_position) != m || any(!is.finite(genetic_position))) {
    .gsim_stop("genetic_position must contain one finite value per variant.")
  }
  if (length(mutation_age) != m || any(!is.finite(mutation_age)) ||
      any(mutation_age < 0)) {
    .gsim_stop(
      "mutation_age must contain one finite nonnegative value per variant."
    )
  }

  chromosome <- as.character(chromosome)
  if (length(chromosome) != m || anyNA(chromosome) || any(!nzchar(chromosome))) {
    .gsim_stop("chromosome must contain one nonmissing label per variant.")
  }
  runs <- rle(chromosome)
  if (anyDuplicated(runs$values)) {
    .gsim_stop("Each chromosome must occupy one contiguous variant block.")
  }
  chromosome_block <- rep.int(seq_along(runs$lengths), runs$lengths)
  block_end <- cumsum(runs$lengths)
  block_start <- c(1L, head(block_end, -1L) + 1L)
  for (i in seq_along(block_start)) {
    idx <- block_start[[i]]:block_end[[i]]
    if (length(idx) > 1L && any(diff(genetic_position[idx]) < 0)) {
      .gsim_stop("genetic_position must be nondecreasing within each chromosome.")
    }
  }

  n <- .gsim_hapnest_integer_scalar(n, "n", 1)
  individual_offset <- .gsim_hapnest_integer_scalar(
    individual_offset, "individual_offset", 0
  )
  if (n > .Machine$integer.max - individual_offset) {
    .gsim_stop("n + individual_offset exceeds the supported range.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) || seed > 2^53 - 1) {
    .gsim_stop("seed must be an integer-valued scalar in [0, 2^53 - 1].")
  }
  if (!is.logical(return_genotypes) || length(return_genotypes) != 1L ||
      is.na(return_genotypes)) {
    .gsim_stop("return_genotypes must be TRUE or FALSE.")
  }

  ans <- .Call(
    C_gsim_hapnest_founders,
    reference_haplotypes_h1,
    reference_haplotypes_h2,
    as.integer(donor_codes),
    weights,
    N,
    Ne,
    rho,
    as.integer(chromosome_block),
    genetic_position,
    mutation_age,
    n,
    as.double(seed),
    return_genotypes,
    individual_offset
  )

  ids <- paste0("syn", individual_offset + seq_len(n))
  variants <- colnames(reference_haplotypes_h1)
  dimnames(ans$h1) <- list(ids, variants)
  dimnames(ans$h2) <- list(ids, variants)
  if (!is.null(ans$genotypes)) dimnames(ans$genotypes) <- list(ids, variants)

  segment <- ans$segments
  segment$chromosome <- runs$values[segment$chromosome_block]
  segment$donor_population <- active[segment$donor_population_code]
  segment$chromosome_block <- NULL
  segment$donor_population_code <- NULL
  segment <- segment[c(
    "individual", "phase", "haplotype", "chromosome", "start", "end",
    "donor_individual", "donor_population", "coalescent_age",
    "sampled_length", "copied_genetic_span", "copied_alternative",
    "retained_alternative"
  )]
  class(segment) <- "data.frame"
  attr(segment, "row.names") <- .set_row_names(length(segment$phase))
  ans$segments <- segment
  ans$settings <- list(
    model = "HAPNEST founder core",
    hapnest_revision = "ba52da1a63cf609306ea92540b3d130fa1efd213",
    donor_phase = donor_phase,
    reference_layout = "paired H1/H2 rows are reference individuals",
    rng = "SplitMix64 per (seed, global haplotype, chromosome block)",
    seed = as.double(seed),
    individual_offset = individual_offset,
    populations = active,
    ancestry_weights = stats::setNames(weights, active),
    N = stats::setNames(N, active),
    Ne = stats::setNames(Ne, active),
    rho = stats::setNames(rho, active),
    chromosomes = runs$values
  )
  class(ans) <- c("gsim_hapnest_founders", "list")
  ans
}

.gsim_hapnest_scales <- function(N, Ne, T, rho) {
  values <- c(N = N, Ne = Ne, T = T, rho = rho)
  if (any(!is.finite(values)) || any(values <= 0)) {
    .gsim_stop("N, Ne, T, and rho must be finite and strictly positive.")
  }
  c(gamma_scale = Ne / N, exponential_scale = 1 / (2 * T * rho))
}

.gsim_hapnest_segment_endpoint <- function(genetic_position, start, last, L) {
  genetic_position <- as.double(genetic_position)
  if (!length(genetic_position) || any(!is.finite(genetic_position)) ||
      any(diff(genetic_position) < 0)) {
    .gsim_stop("genetic_position must be a finite nondecreasing vector.")
  }
  start <- .gsim_hapnest_integer_scalar(start, "start", 1)
  last <- .gsim_hapnest_integer_scalar(last, "last", 1)
  if (start > last || last > length(genetic_position)) {
    .gsim_stop("Require 1 <= start <= last <= length(genetic_position).")
  }
  if (!is.numeric(L) || length(L) != 1L || is.na(L) || L < 0) {
    .gsim_stop("L must be a nonnegative scalar.")
  }
  .Call(C_gsim_hapnest_segment_endpoint, genetic_position, start, last,
        as.double(L))
}

.gsim_hapnest_copy_segment <- function(donor, mutation_age, start, end, T) {
  donor <- as.vector(.gsim_hapnest_raw_matrix(matrix(donor, nrow = 1L), "donor"))
  mutation_age <- as.double(mutation_age)
  if (length(mutation_age) != length(donor) || any(!is.finite(mutation_age)) ||
      any(mutation_age < 0)) {
    .gsim_stop("mutation_age must be finite, nonnegative, and aligned to donor.")
  }
  start <- .gsim_hapnest_integer_scalar(start, "start", 1)
  end <- .gsim_hapnest_integer_scalar(end, "end", 1)
  if (start > end || end > length(donor)) {
    .gsim_stop("Require 1 <= start <= end <= length(donor).")
  }
  if (!is.numeric(T) || length(T) != 1L || is.na(T) || !is.finite(T) || T < 0) {
    .gsim_stop("T must be a finite nonnegative scalar.")
  }
  .Call(C_gsim_hapnest_copy_segment, donor, mutation_age, start, end,
        as.double(T))
}

.gsim_hapnest_pair <- function(h1, h2) {
  h1 <- .gsim_hapnest_raw_matrix(h1, "h1")
  h2 <- .gsim_hapnest_raw_matrix(h2, "h2")
  if (!identical(dim(h1), dim(h2))) .gsim_stop("h1 and h2 must have equal dimensions.")
  ans <- .Call(C_gsim_hapnest_pair, h1, h2)
  dimnames(ans) <- dimnames(h1)
  ans
}
