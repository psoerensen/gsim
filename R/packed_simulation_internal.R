# Private chromosome-local packed backend compiled into gsim.
.gsim_packed_backend <- function() {
  pointer <- .Call(C_gsim_packed_backend)
  attr(pointer, "packed_origin") <- "internalized gbits 0.20.0 (089bf1e)"
  attr(pointer, "packed_abi") <- NA_integer_
  class(pointer) <- "gsim_packed_backend"
  pointer
}

.gsim_packed_tag <- function(pointer, sample_ids = NULL, variant_ids = NULL) {
  attr(pointer, "sample_ids") <- sample_ids
  attr(pointer, "variant_ids") <- variant_ids
  class(pointer) <- "gsim_packed_haplotypes"
  pointer
}

.gsim_packed_pack <- function(backend, values) {
  values <- .gsim_hapnest_raw_matrix(values, "values")
  .gsim_packed_tag(
    .Call(C_gsim_packed_pack, backend, values), rownames(values), colnames(values)
  )
}

.gsim_packed_zero <- function(backend, individuals, markers,
                             sample_ids = NULL, variant_ids = NULL) {
  individuals <- .gsim_hapnest_integer_scalar(individuals, "individuals", 1)
  markers <- .gsim_hapnest_integer_scalar(markers, "markers", 1)
  .gsim_packed_tag(
    .Call(C_gsim_packed_zero, backend, individuals, markers),
    sample_ids, variant_ids
  )
}

.gsim_packed_set_marker <- function(h1, h2, marker, h1_values, h2_values) {
  invisible(.Call(C_gsim_packed_set_marker, h1, h2, as.integer(marker),
                  h1_values, h2_values))
}

.gsim_packed_unpack <- function(haplotypes) {
  out <- .Call(C_gsim_packed_unpack, haplotypes)
  sample_ids <- attr(haplotypes, "sample_ids", exact = TRUE)
  variant_ids <- attr(haplotypes, "variant_ids", exact = TRUE)
  if (!is.null(sample_ids) || !is.null(variant_ids)) {
    dimnames(out) <- list(sample_ids, variant_ids)
  }
  out
}

.gsim_packed_info <- function(haplotypes) {
  .Call(C_gsim_packed_info, haplotypes)
}

.gsim_packed_close <- function(haplotypes) {
  invisible(.Call(C_gsim_packed_close, haplotypes))
}

.gsim_packed_word <- function(haplotypes, marker, word) {
  marker <- .gsim_hapnest_integer_scalar(marker, "marker", 1)
  word <- .gsim_hapnest_integer_scalar(word, "word", 1)
  .Call(C_gsim_packed_word, haplotypes, marker - 1L, word - 1L)
}

.gsim_packed_copy_interval <- function(destination, destination_individual,
                                      source, source_individual, first, last) {
  values <- c(
    destination_individual = .gsim_hapnest_integer_scalar(
      destination_individual, "destination_individual", 1
    ),
    source_individual = .gsim_hapnest_integer_scalar(
      source_individual, "source_individual", 1
    ),
    first = .gsim_hapnest_integer_scalar(first, "first", 1),
    last = .gsim_hapnest_integer_scalar(last, "last", 1)
  )
  invisible(.Call(
    C_gsim_packed_copy_interval, destination, values[[1L]] - 1L,
    source, values[[2L]] - 1L, values[[3L]] - 1L, values[[4L]] - 1L
  ))
}

.gsim_packed_copy_filtered <- function(destination, destination_individual,
                                      source, source_individual, first, last,
                                      coalescent_age, mutation_age) {
  if (!is.numeric(coalescent_age) || length(coalescent_age) != 1L ||
      !is.finite(coalescent_age) || coalescent_age < 0) {
    .gsim_stop("coalescent_age must be one finite nonnegative value.")
  }
  mutation_age <- as.double(mutation_age)
  invisible(.Call(
    C_gsim_packed_copy_filtered, destination,
    as.integer(destination_individual - 1L), source,
    as.integer(source_individual - 1L), as.integer(first - 1L),
    as.integer(last - 1L), as.double(coalescent_age), mutation_age
  ))
}

.gsim_packed_copy_filtered_counts <- function(
  destination, destination_individual, source, source_individual, first, last,
  coalescent_age, mutation_age
) {
  if (!is.numeric(coalescent_age) || length(coalescent_age) != 1L ||
      !is.finite(coalescent_age) || coalescent_age < 0) {
    .gsim_stop("coalescent_age must be one finite nonnegative value.")
  }
  .Call(
    C_gsim_packed_copy_filtered_counts, destination,
    as.integer(destination_individual - 1L), source,
    as.integer(source_individual - 1L), as.integer(first - 1L),
    as.integer(last - 1L), as.double(coalescent_age), as.double(mutation_age)
  )
}

.gsim_packed_make_gamete <- function(destination, destination_individual,
                                    parent_h1, parent_h2, parent_individual,
                                    starting_haplotype, boundaries) {
  invisible(.Call(
    C_gsim_packed_make_gamete, destination,
    as.integer(destination_individual - 1L), parent_h1, parent_h2,
    as.integer(parent_individual - 1L), as.integer(starting_haplotype),
    as.integer(boundaries)
  ))
}

.gsim_packed_decode_genotypes <- function(h1, h2) {
  out <- .Call(C_gsim_packed_decode_genotypes, h1, h2)
  dimnames(out) <- list(
    attr(h1, "sample_ids", exact = TRUE),
    attr(h1, "variant_ids", exact = TRUE)
  )
  out
}

.gsim_hapnest_packed_reference_inputs <- function(
  backend,
  reference_h1,
  reference_h2,
  donor_population,
  ancestry_weights,
  N,
  Ne,
  rho,
  genetic_position,
  mutation_age,
  n,
  seed,
  chromosome,
  donor_phase = "hapnest",
  return_genotypes = FALSE,
  return_segments = TRUE,
  individual_offset = 0L
) {
  if (!inherits(backend, "gsim_packed_backend") ||
      !inherits(reference_h1, "gsim_packed_haplotypes") ||
      !inherits(reference_h2, "gsim_packed_haplotypes")) {
    .gsim_stop("backend and packed reference H1/H2 handles are required.")
  }
  info_h1 <- .gsim_packed_info(reference_h1)
  info_h2 <- .gsim_packed_info(reference_h2)
  if (!identical(unname(info_h1[1:3]), unname(info_h2[1:3]))) {
    .gsim_stop("Packed reference H1 and H2 dimensions and layout must match.")
  }
  donor_count <- as.integer(info_h1[[1L]])
  marker_count <- as.integer(info_h1[[2L]])
  donor_ids <- attr(reference_h1, "sample_ids", exact = TRUE)
  variant_ids <- attr(reference_h1, "variant_ids", exact = TRUE)
  if (is.null(donor_ids) || length(donor_ids) != donor_count || anyNA(donor_ids) ||
      any(!nzchar(donor_ids)) || anyDuplicated(donor_ids) ||
      !identical(donor_ids, attr(reference_h2, "sample_ids", exact = TRUE))) {
    .gsim_stop("Packed reference H1/H2 require identical unique donor IDs.")
  }
  if (is.null(variant_ids) || length(variant_ids) != marker_count ||
      anyNA(variant_ids) || any(!nzchar(variant_ids)) || anyDuplicated(variant_ids) ||
      !identical(variant_ids, attr(reference_h2, "variant_ids", exact = TRUE))) {
    .gsim_stop("Packed reference H1/H2 require identical unique variant IDs.")
  }
  donor_ids <- enc2utf8(as.character(donor_ids))
  variant_ids <- enc2utf8(as.character(variant_ids))
  if (!is.character(donor_phase) || length(donor_phase) != 1L ||
      is.na(donor_phase) || donor_phase != "hapnest") {
    .gsim_stop("donor_phase must be 'hapnest'; pooled sampling is not implemented.")
  }
  donor_population <- as.character(donor_population)
  if (!is.null(names(donor_population))) {
    population_ids <- enc2utf8(names(donor_population))
    if (anyNA(population_ids) || any(!nzchar(population_ids)) ||
        anyDuplicated(population_ids) || !setequal(population_ids, donor_ids) ||
        length(population_ids) != donor_count) {
      .gsim_stop("Named donor_population IDs must exactly match packed donor IDs.")
    }
    donor_population <- donor_population[match(donor_ids, population_ids)]
  }
  if (length(donor_population) != donor_count || anyNA(donor_population) ||
      any(!nzchar(donor_population))) {
    .gsim_stop("donor_population must align one nonempty label per packed donor.")
  }
  if (!is.numeric(ancestry_weights) || is.null(names(ancestry_weights)) ||
      anyNA(names(ancestry_weights)) || any(!nzchar(names(ancestry_weights))) ||
      anyDuplicated(names(ancestry_weights)) ||
      any(!is.finite(ancestry_weights)) || any(ancestry_weights < 0)) {
    .gsim_stop("ancestry_weights must be uniquely named, finite, and nonnegative.")
  }
  active <- names(ancestry_weights)[ancestry_weights > 0]
  if (!length(active) || !is.finite(sum(ancestry_weights[active]))) {
    .gsim_stop("ancestry_weights must have a positive finite sum.")
  }
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
  observed_N <- vapply(seq_along(active), function(i) sum(donor_codes == i),
                       integer(1L))
  if (!identical(N, as.double(observed_N))) {
    .gsim_stop("N must equal packed reference individuals in each active population.")
  }
  genetic_position <- as.double(genetic_position)
  mutation_age <- as.double(mutation_age)
  if (length(genetic_position) != marker_count ||
      any(!is.finite(genetic_position)) || any(diff(genetic_position) < 0)) {
    .gsim_stop("genetic_position must be finite, nondecreasing, and marker-aligned.")
  }
  if (length(mutation_age) != marker_count || any(!is.finite(mutation_age)) ||
      any(mutation_age < 0)) {
    .gsim_stop("mutation_age must be finite, nonnegative, and marker-aligned.")
  }
  chromosome <- as.character(chromosome)
  if (length(chromosome) == 1L) chromosome <- rep.int(chromosome, marker_count)
  if (length(chromosome) != marker_count || anyNA(chromosome) ||
      any(!nzchar(chromosome)) || length(unique(chromosome)) != 1L) {
    .gsim_stop("The packed founder call requires exactly one chromosome label.")
  }
  n <- .gsim_hapnest_integer_scalar(n, "n", 1)
  individual_offset <- .gsim_hapnest_integer_scalar(
    individual_offset, "individual_offset", 0)
  if (n > .Machine$integer.max - individual_offset) {
    .gsim_stop("n + individual_offset exceeds the supported range.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) || seed > 2^53 - 1) {
    .gsim_stop("seed must be an integer-valued scalar in [0, 2^53 - 1].")
  }
  return_genotypes <- .gsim_meiosis_flag(return_genotypes, "return_genotypes")
  return_segments <- .gsim_meiosis_flag(return_segments, "return_segments")
  list(
    info_h1 = info_h1, info_h2 = info_h2, donor_count = donor_count,
    marker_count = marker_count, donor_ids = donor_ids, variant_ids = variant_ids,
    donor_codes = as.integer(donor_codes), active = active, weights = weights,
    N = N, Ne = Ne, rho = rho, genetic_position = genetic_position,
    mutation_age = mutation_age, n = n, seed = as.double(seed),
    chromosome = chromosome, donor_phase = donor_phase,
    return_genotypes = return_genotypes, return_segments = return_segments,
    individual_offset = individual_offset
  )
}

.gsim_hapnest_founders_packed_reference_chromosome <- function(
  backend, reference_h1, reference_h2, donor_population, ancestry_weights,
  N, Ne, rho, genetic_position, mutation_age, n, seed, chromosome,
  donor_phase = "hapnest", return_genotypes = FALSE, return_segments = TRUE,
  individual_offset = 0L
) {
  input <- .gsim_hapnest_packed_reference_inputs(
    backend, reference_h1, reference_h2, donor_population, ancestry_weights,
    N, Ne, rho, genetic_position, mutation_age, n, seed, chromosome,
    donor_phase, return_genotypes, return_segments, individual_offset)
  segment <- .Call(
    C_gsim_hapnest_plan, input$donor_count, input$marker_count,
    input$donor_codes, input$weights, input$N, input$Ne, input$rho,
    rep.int(1L, input$marker_count), input$chromosome[[1L]],
    input$genetic_position, input$n, input$seed, input$individual_offset)
  segment$chromosome <- rep.int(input$chromosome[[1L]], length(segment$phase))
  segment$donor_population <- input$active[segment$donor_population_code]
  segment$chromosome_block <- NULL
  segment$donor_population_code <- NULL
  segment <- segment[c(
    "individual", "phase", "haplotype", "chromosome", "start", "end",
    "donor_individual", "donor_population", "coalescent_age",
    "sampled_length", "copied_genetic_span", "copied_alternative",
    "retained_alternative")]
  class(segment) <- "data.frame"
  attr(segment, "row.names") <- .set_row_names(length(segment$phase))
  ids <- paste0("syn", input$individual_offset + seq_len(input$n))
  h1 <- .gsim_packed_zero(backend, input$n, input$marker_count,
                         ids, input$variant_ids)
  h2 <- .gsim_packed_zero(backend, input$n, input$marker_count,
                         ids, input$variant_ids)
  materialized <- FALSE
  on.exit({
    if (!materialized) {
      try(.gsim_packed_close(h1), silent = TRUE)
      try(.gsim_packed_close(h2), silent = TRUE)
    }
  }, add = TRUE)
  for (i in seq_len(nrow(segment))) {
    record <- segment[i, ]
    destination <- if (record$phase == 1L) h1 else h2
    source <- if (record$phase == 1L) reference_h1 else reference_h2
    if (input$return_segments) {
      counts <- .gsim_packed_copy_filtered_counts(
        destination, record$individual - input$individual_offset, source,
        record$donor_individual, record$start, record$end,
        record$coalescent_age, input$mutation_age)
      segment$copied_alternative[[i]] <- counts[[1L]]
      segment$retained_alternative[[i]] <- counts[[2L]]
    } else {
      .gsim_packed_copy_filtered(
        destination, record$individual - input$individual_offset, source,
        record$donor_individual, record$start, record$end,
        record$coalescent_age, input$mutation_age)
    }
  }
  genotypes <- if (input$return_genotypes) {
    .gsim_packed_decode_genotypes(h1, h2)
  } else NULL
  settings <- list(
    model = "HAPNEST founder core",
    hapnest_revision = "ba52da1a63cf609306ea92540b3d130fa1efd213",
    donor_phase = "hapnest",
    reference_layout = "paired packed H1/H2 rows are reference individuals",
    reference_source = "caller-owned gsim packed handles",
    chromosome_identity = "exact UTF-8 label bytes; no normalization",
    chromosome_hash = "FNV-1a 64-bit",
    rng = "SplitMix64 per (seed, global haplotype, chromosome identity)",
    seed = input$seed, individual_offset = input$individual_offset,
    populations = input$active,
    ancestry_weights = stats::setNames(input$weights, input$active),
    N = stats::setNames(input$N, input$active),
    Ne = stats::setNames(input$Ne, input$active),
    rho = stats::setNames(input$rho, input$active),
    chromosomes = input$chromosome[[1L]],
    storage = "gsim marker-major one-bit phased haplotypes",
    implementation_origin = attr(backend, "packed_origin", exact = TRUE),
    decoded_genotypes = input$return_genotypes)
  reference_bytes <- unname(input$info_h1[[4L]] + input$info_h2[[4L]])
  generated_bytes <- unname(.gsim_packed_info(h1)[[4L]] +
                              .gsim_packed_info(h2)[[4L]])
  out <- list(
    h1 = h1, h2 = h2, genotypes = genotypes,
    segments = if (input$return_segments) segment else NULL,
    sample_ids = ids, variant_ids = input$variant_ids,
    chromosome = input$chromosome[[1L]], settings = settings,
    memory = list(
      reference_raw_bytes_avoided = 2 * input$donor_count * input$marker_count,
      reference_packed_bytes = reference_bytes,
      generated_raw_bytes_avoided = 2 * input$n * input$marker_count,
      generated_packed_bytes = generated_bytes,
      event_record_bytes = as.numeric(object.size(segment)),
      peak_biological_payload_bytes = reference_bytes + generated_bytes,
      chromosome_temporary_bytes = 0,
      decoded_genotype_bytes = if (input$return_genotypes) length(genotypes) else 0)
  )
  class(out) <- c("gsim_hapnest_founders_packed", "list")
  materialized <- TRUE
  out
}

.gsim_hapnest_founders_packed_chromosome <- function(
  backend, reference_haplotypes_h1, reference_haplotypes_h2,
  donor_population, ancestry_weights, N, Ne, rho, genetic_position,
  mutation_age, n, seed, chromosome, donor_phase = "hapnest",
  return_genotypes = FALSE, return_segments = TRUE, individual_offset = 0L
) {
  raw_h1 <- .gsim_hapnest_raw_matrix(
    reference_haplotypes_h1, "reference_haplotypes_h1")
  raw_h2 <- .gsim_hapnest_raw_matrix(
    reference_haplotypes_h2, "reference_haplotypes_h2")
  if (!identical(dim(raw_h1), dim(raw_h2)) ||
      !identical(dimnames(raw_h1), dimnames(raw_h2))) {
    .gsim_stop("Reference H1/H2 dimensions and names must be identical.")
  }
  reference_h1 <- .gsim_packed_pack(backend, raw_h1)
  reference_h2 <- .gsim_packed_pack(backend, raw_h2)
  on.exit({
    try(.gsim_packed_close(reference_h1), silent = TRUE)
    try(.gsim_packed_close(reference_h2), silent = TRUE)
  }, add = TRUE)
  out <- .gsim_hapnest_founders_packed_reference_chromosome(
    backend, reference_h1, reference_h2, donor_population, ancestry_weights,
    N, Ne, rho, genetic_position, mutation_age, n, seed, chromosome,
    donor_phase, return_genotypes, return_segments, individual_offset)
  out$memory$reference_raw_bytes <- 2 * length(raw_h1)
  out
}

.gsim_pedigree_genotypes_packed_chromosome <- function(
  backend,
  pedigree,
  founder_haplotypes,
  chromosome,
  genetic_position,
  seed,
  return_haplotypes = TRUE,
  return_genotypes = FALSE,
  return_crossovers = TRUE,
  batch_size = NULL
) {
  if (!inherits(pedigree, "gsim_pedigree") ||
      !is.data.frame(pedigree$pedigree) || is.null(pedigree$canonical_order)) {
    .gsim_stop("pedigree must be a gsim_pedigree object.")
  }
  canonical <- as.character(pedigree$canonical_order)
  tab <- pedigree$pedigree
  required <- c("animal", "sire", "dam")
  if (any(!required %in% names(tab)) || !length(canonical) ||
      anyNA(canonical) || any(!nzchar(canonical)) || anyDuplicated(canonical) ||
      nrow(tab) != length(canonical) || anyNA(tab$animal) ||
      any(!nzchar(as.character(tab$animal))) || anyDuplicated(tab$animal) ||
      !setequal(as.character(tab$animal), canonical)) {
    .gsim_stop("pedigree canonical and external animal mappings are invalid.")
  }
  tab <- tab[match(canonical, as.character(tab$animal)), , drop = FALSE]
  rownames(tab) <- NULL
  tab$animal <- as.character(tab$animal)
  tab$sire <- as.character(tab$sire)
  tab$dam <- as.character(tab$dam)
  sire_position <- match(tab$sire, canonical)
  dam_position <- match(tab$dam, canonical)
  founder <- is.na(tab$sire) & is.na(tab$dam)
  partial <- xor(is.na(tab$sire), is.na(tab$dam))
  if (any(partial)) {
    i <- which(partial)[[1L]]
    missing_side <- if (is.na(tab$sire[[i]])) "sire" else "dam"
    .gsim_stop("Animal '", tab$animal[[i]],
               "' has one known parent; missing ", missing_side,
               " is not supported.")
  }
  for (i in which(!founder)) {
    if (is.na(sire_position[[i]]) || is.na(dam_position[[i]])) {
      .gsim_stop("Recorded parent of animal '", tab$animal[[i]],
                 "' is absent from canonical_order.")
    }
    if (sire_position[[i]] >= i || dam_position[[i]] >= i) {
      .gsim_stop("Canonical parent-before-offspring contract is violated for animal '",
                 tab$animal[[i]], "'.")
    }
  }

  if (!is.list(founder_haplotypes) ||
      !all(c("h1", "h2") %in% names(founder_haplotypes))) {
    .gsim_stop("founder_haplotypes must be a named list containing h1 and h2.")
  }
  packed_founders <- inherits(founder_haplotypes$h1, "gsim_packed_haplotypes") &&
    inherits(founder_haplotypes$h2, "gsim_packed_haplotypes")
  if (xor(inherits(founder_haplotypes$h1, "gsim_packed_haplotypes"),
          inherits(founder_haplotypes$h2, "gsim_packed_haplotypes"))) {
    .gsim_stop("Founder H1 and H2 must both be packed handles or both be matrices.")
  }
  if (packed_founders) {
    source_h1 <- founder_haplotypes$h1
    source_h2 <- founder_haplotypes$h2
    info_h1 <- .gsim_packed_info(source_h1)
    info_h2 <- .gsim_packed_info(source_h2)
    if (!identical(unname(info_h1[1:3]), unname(info_h2[1:3]))) {
      .gsim_stop("Founder packed H1 and H2 dimensions and layout must match.")
    }
    supplied_ids <- attr(source_h1, "sample_ids", exact = TRUE)
    variant_ids <- attr(source_h1, "variant_ids", exact = TRUE)
    if (is.null(supplied_ids) || anyNA(supplied_ids) ||
        any(!nzchar(supplied_ids)) || anyDuplicated(supplied_ids) ||
        !identical(supplied_ids,
                   attr(source_h2, "sample_ids", exact = TRUE))) {
      .gsim_stop("Founder packed H1/H2 require identical unique animal IDs.")
    }
    if (is.null(variant_ids) || anyNA(variant_ids) ||
        any(!nzchar(variant_ids)) || anyDuplicated(variant_ids) ||
        !identical(variant_ids,
                   attr(source_h2, "variant_ids", exact = TRUE))) {
      .gsim_stop("Founder packed H1/H2 require identical unique variant IDs.")
    }
    founder_row_count <- as.integer(info_h1[[1L]])
    marker_count <- as.integer(info_h1[[2L]])
    if (length(supplied_ids) != founder_row_count ||
        length(variant_ids) != marker_count) {
      .gsim_stop("Founder packed IDs must exactly match handle dimensions.")
    }
    founder_raw_bytes <- 2 * founder_row_count * marker_count
  } else {
    founder_h1 <- .gsim_hapnest_raw_matrix(
      founder_haplotypes$h1, "founder_haplotypes$h1")
    founder_h2 <- .gsim_hapnest_raw_matrix(
      founder_haplotypes$h2, "founder_haplotypes$h2")
    if (!identical(dim(founder_h1), dim(founder_h2))) {
      .gsim_stop("Founder H1 and H2 matrices must have identical dimensions.")
    }
    if (is.null(rownames(founder_h1)) || is.null(rownames(founder_h2)) ||
        anyNA(rownames(founder_h1)) || any(!nzchar(rownames(founder_h1))) ||
        anyDuplicated(rownames(founder_h1)) ||
        !identical(rownames(founder_h1), rownames(founder_h2))) {
      .gsim_stop("Founder H1/H2 row names must be identical unique animal IDs.")
    }
    if (is.null(colnames(founder_h1)) || is.null(colnames(founder_h2)) ||
        anyNA(colnames(founder_h1)) || any(!nzchar(colnames(founder_h1))) ||
        anyDuplicated(colnames(founder_h1)) ||
        !identical(colnames(founder_h1), colnames(founder_h2))) {
      .gsim_stop("Founder H1/H2 column names must be identical unique variant IDs.")
    }
    supplied_ids <- rownames(founder_h1)
    variant_ids <- colnames(founder_h1)
    founder_row_count <- nrow(founder_h1)
    marker_count <- ncol(founder_h1)
    founder_raw_bytes <- 2 * length(founder_h1)
    source_h1 <- .gsim_packed_pack(backend, founder_h1)
    source_h2 <- .gsim_packed_pack(backend, founder_h2)
    on.exit({
      try(.gsim_packed_close(source_h1), silent = TRUE)
      try(.gsim_packed_close(source_h2), silent = TRUE)
    }, add = TRUE)
  }
  founder_ids <- canonical[founder]
  missing_founders <- setdiff(founder_ids, supplied_ids)
  extra_founders <- setdiff(supplied_ids, founder_ids)
  if (length(missing_founders) || length(extra_founders) ||
      founder_row_count != length(founder_ids)) {
    detail <- c(
      if (length(missing_founders)) paste0("missing: ", paste(missing_founders, collapse = ", ")),
      if (length(extra_founders)) paste0("extra: ", paste(extra_founders, collapse = ", "))
    )
    .gsim_stop("Founder haplotype IDs must exactly match pedigree founders",
               if (length(detail)) paste0(" (", paste(detail, collapse = "; "), ")") else "",
               ".")
  }
  chromosome <- as.character(chromosome)
  if (length(chromosome) != marker_count || anyNA(chromosome) ||
      any(!nzchar(chromosome)) || length(unique(chromosome)) != 1L) {
    .gsim_stop("The packed pedigree call requires one chromosome label per variant and exactly one chromosome.")
  }
  genetic_position <- as.double(genetic_position)
  if (length(genetic_position) != marker_count ||
      any(!is.finite(genetic_position)) || any(diff(genetic_position) < 0)) {
    .gsim_stop("genetic_position must be finite, nondecreasing, and aligned.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) || seed > 2^53 - 1) {
    .gsim_stop("seed must be an integer-valued scalar in [0, 2^53 - 1].")
  }
  return_haplotypes <- .gsim_meiosis_flag(return_haplotypes, "return_haplotypes")
  return_genotypes <- .gsim_meiosis_flag(return_genotypes, "return_genotypes")
  return_crossovers <- .gsim_meiosis_flag(return_crossovers, "return_crossovers")
  if (is.null(batch_size)) {
    batch_size <- length(canonical)
  } else {
    batch_size <- .gsim_hapnest_integer_scalar(batch_size, "batch_size", 1)
  }

  animal_count <- length(canonical)
  h1 <- .gsim_packed_zero(backend, animal_count, marker_count,
                         canonical, variant_ids)
  h2 <- .gsim_packed_zero(backend, animal_count, marker_count,
                         canonical, variant_ids)
  founder_order <- match(founder_ids, supplied_ids)
  founder_positions <- match(founder_ids, canonical)
  for (i in seq_along(founder_ids)) {
    .gsim_packed_copy_interval(h1, founder_positions[[i]], source_h1,
                              founder_order[[i]], 1L, marker_count)
    .gsim_packed_copy_interval(h2, founder_positions[[i]], source_h2,
                              founder_order[[i]], 1L, marker_count)
  }

  meiosis_rows <- list()
  crossover_rows <- list()
  meiosis_number <- 0L
  crossover_number <- 0L
  batch_starts <- seq.int(1L, animal_count, by = batch_size)
  for (batch_start in batch_starts) {
    batch_end <- min(animal_count, batch_start + batch_size - 1L)
    for (i in seq.int(batch_start, batch_end)) {
      if (founder[[i]]) next
      for (side in 1:2) {
        parent_position <- if (side == 1L) sire_position[[i]] else dam_position[[i]]
        parent_id <- canonical[[parent_position]]
        plan <- .Call(
          C_gsim_meiosis_plan, genetic_position, as.double(seed), canonical[[i]],
          chromosome[[1L]], as.integer(side)
        )
        destination <- if (side == 1L) h1 else h2
        .gsim_packed_make_gamete(
          destination, i, h1, h2, parent_position,
          plan$starting_haplotype, plan$boundaries
        )
        if (return_crossovers) {
          side_name <- if (side == 1L) "paternal" else "maternal"
          child_phase <- if (side == 1L) "H1" else "H2"
          meiosis_number <- meiosis_number + 1L
          meiosis_rows[[meiosis_number]] <- data.frame(
            animal = canonical[[i]], parent = parent_id,
            parental_side = side_name, child_phase = child_phase,
            chromosome = chromosome[[1L]],
            starting_haplotype = plan$starting_haplotype,
            crossover_count = length(plan$crossovers),
            stringsAsFactors = FALSE
          )
          if (length(plan$crossovers)) {
            crossover_number <- crossover_number + 1L
            crossover_rows[[crossover_number]] <- data.frame(
              animal = rep.int(canonical[[i]], length(plan$crossovers)),
              parent = rep.int(parent_id, length(plan$crossovers)),
              parental_side = rep.int(side_name, length(plan$crossovers)),
              child_phase = rep.int(child_phase, length(plan$crossovers)),
              chromosome = rep.int(chromosome[[1L]], length(plan$crossovers)),
              crossover_index = seq_along(plan$crossovers),
              genetic_position = plan$crossovers,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  genotypes <- if (return_genotypes) .gsim_packed_decode_genotypes(h1, h2) else NULL
  crossover_audit <- NULL
  if (return_crossovers) {
    crossover_audit <- list(
      meioses = if (length(meiosis_rows)) do.call(rbind, meiosis_rows) else
        .gsim_meiosis_empty_meioses(),
      crossovers = if (length(crossover_rows)) do.call(rbind, crossover_rows) else
        .gsim_meiosis_empty_crossovers()
    )
    rownames(crossover_audit$meioses) <- NULL
    rownames(crossover_audit$crossovers) <- NULL
  }
  pedigree_alignment <- data.frame(
    animal = canonical, canonical_index = seq_along(canonical),
    external_index = match(canonical, as.character(pedigree$pedigree$animal)),
    sire = tab$sire, dam = tab$dam, founder = founder,
    stringsAsFactors = FALSE
  )
  variant_map <- data.frame(
    variant = variant_ids, variant_index = seq_along(variant_ids),
    chromosome = chromosome, genetic_position = genetic_position,
    stringsAsFactors = FALSE
  )
  meiosis_event_bytes <- if (return_crossovers) {
    as.numeric(object.size(crossover_audit))
  } else {
    0
  }
  out <- list(
    h1 = if (return_haplotypes) h1 else NULL,
    h2 = if (return_haplotypes) h2 else NULL,
    genotypes = genotypes, crossover_audit = crossover_audit,
    sample_ids = canonical, variant_ids = variant_ids,
    variant_map = variant_map,
    chromosome_blocks = data.frame(
      chromosome = chromosome[[1L]], start_variant = 1L,
      end_variant = marker_count,
      length_morgans = genetic_position[[marker_count]] - genetic_position[[1L]],
      stringsAsFactors = FALSE
    ),
    pedigree_alignment = pedigree_alignment,
    settings = list(
      model = "marker-level no-interference Mendelian meiosis",
      phase = c(H1 = "paternal gamete", H2 = "maternal gamete"),
      founder_phase = "supplied H1/H2 labels retained",
      genetic_position_unit = "Morgan",
      crossover_count = "Poisson(last_position - first_position)",
      crossover_location = "sorted conditional Uniform(first_position, last_position)",
      crossover_boundary = "crossovers at x switch before assigning markers at x",
      rng = "SplitMix64 with stable FNV-1a UTF-8 identifier hashing",
      stream_key = "seed, child animal ID, chromosome label, parental side",
      seed = as.double(seed), global_rng = "not used",
      batch_policy = "operational batch boundaries do not enter RNG streams",
      unrelated_animal_policy = "stable IDs leave existing streams unchanged",
      storage = "gsim marker-major one-bit phased haplotypes",
      founder_source = if (packed_founders) "caller-owned packed handles" else
        "caller-supplied byte matrices packed internally",
      implementation_origin = attr(backend, "packed_origin", exact = TRUE),
      decoded_genotypes = return_genotypes
    ),
    memory = list(
      founder_raw_bytes = founder_raw_bytes,
      founder_packed_bytes = unname(.gsim_packed_info(source_h1)[[4L]] +
                                      .gsim_packed_info(source_h2)[[4L]]),
      generated_raw_bytes = 2 * animal_count * marker_count,
      generated_packed_bytes = unname(.gsim_packed_info(h1)[[4L]] +
                                        .gsim_packed_info(h2)[[4L]]),
      event_record_bytes = meiosis_event_bytes,
      chromosome_temporary_bytes = 8 * max(
        1L, if (return_crossovers && nrow(crossover_audit$meioses))
          max(crossover_audit$meioses$crossover_count) else 0L
      ),
      decoded_genotype_bytes = if (return_genotypes) length(genotypes) else 0
    )
  )
  class(out) <- c("gsim_pedigree_genotypes_packed", "list")
  out
}
