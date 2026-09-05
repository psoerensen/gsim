# Experimental chromosome-local gbits integration.  This interface is
# intentionally unexported; statistical policy and event generation remain in
# gsim, while the dynamically resolved gbits ABI owns packed storage.

.gsim_gbits_symbols <- c(
  "gbits_abi_version", "gbits_library_version", "gbits_last_error",
  "gbits_phased_haplotype_create_zero",
  "gbits_phased_haplotype_create_from_values",
  "gbits_phased_haplotype_close",
  "gbits_phased_haplotype_individual_count",
  "gbits_phased_haplotype_marker_count",
  "gbits_phased_haplotype_words_per_marker",
  "gbits_phased_haplotype_storage_bytes",
  "gbits_phased_haplotype_word",
  "gbits_phased_haplotype_allele",
  "gbits_phased_haplotype_set_allele",
  "gbits_phased_haplotype_unpack",
  "gbits_phased_haplotype_copy_interval",
  "gbits_phased_haplotype_copy_filtered_segment",
  "gbits_phased_haplotype_make_gamete",
  "gbits_phased_haplotype_decode_genotypes",
  "gbits_bed_open", "gbits_bed_close", "gbits_bed_read_variant",
  "gbits_bed_sink_create", "gbits_bed_sink_append_phased",
  "gbits_bed_sink_finalize", "gbits_bed_sink_get_info",
  "gbits_bed_sink_close",
  "gbits_hap_sink_create", "gbits_hap_sink_append_phased",
  "gbits_hap_sink_finalize", "gbits_hap_sink_get_info",
  "gbits_hap_sink_close", "gbits_hap_reader_open",
  "gbits_hap_reader_close", "gbits_hap_reader_dimensions",
  "gbits_hap_reader_chromosome_info",
  "gbits_hap_reader_load_chromosome"
)

.gsim_gbits_backend <- function(
  library = Sys.getenv("GSIM_GBITS_LIBRARY", unset = "")
) {
  if (!is.character(library) || length(library) != 1L || is.na(library) ||
      !nzchar(library)) {
    .gsim_stop(
      "A gbits shared-library path is required via library or GSIM_GBITS_LIBRARY."
    )
  }
  path <- normalizePath(library, winslash = "/", mustWork = TRUE)
  dll <- tryCatch(
    dyn.load(path, local = TRUE, now = TRUE),
    error = function(error) .gsim_stop("Unable to load gbits: ", conditionMessage(error))
  )
  addresses <- tryCatch(
    stats::setNames(lapply(.gsim_gbits_symbols, function(name) {
      getNativeSymbolInfo(name, PACKAGE = dll)$address
    }), .gsim_gbits_symbols),
    error = function(error) {
      .gsim_stop("Incompatible gbits library: ", conditionMessage(error))
    }
  )
  pointer <- .Call(C_gsim_gbits_backend, addresses)
  attr(pointer, "library") <- path
  attr(pointer, "dll") <- dll
  attr(pointer, "gbits_abi") <- 4L
  class(pointer) <- "gsim_gbits_backend"
  pointer
}

.gsim_gbits_tag <- function(pointer, sample_ids = NULL, variant_ids = NULL) {
  attr(pointer, "sample_ids") <- sample_ids
  attr(pointer, "variant_ids") <- variant_ids
  class(pointer) <- "gsim_gbits_haplotypes"
  pointer
}

.gsim_gbits_pack <- function(backend, values) {
  values <- .gsim_hapnest_raw_matrix(values, "values")
  .gsim_gbits_tag(
    .Call(C_gsim_gbits_pack, backend, values), rownames(values), colnames(values)
  )
}

.gsim_gbits_zero <- function(backend, individuals, markers,
                             sample_ids = NULL, variant_ids = NULL) {
  individuals <- .gsim_hapnest_integer_scalar(individuals, "individuals", 1)
  markers <- .gsim_hapnest_integer_scalar(markers, "markers", 1)
  .gsim_gbits_tag(
    .Call(C_gsim_gbits_zero, backend, individuals, markers),
    sample_ids, variant_ids
  )
}

.gsim_gbits_unpack <- function(haplotypes) {
  out <- .Call(C_gsim_gbits_unpack, haplotypes)
  sample_ids <- attr(haplotypes, "sample_ids", exact = TRUE)
  variant_ids <- attr(haplotypes, "variant_ids", exact = TRUE)
  if (!is.null(sample_ids) || !is.null(variant_ids)) {
    dimnames(out) <- list(sample_ids, variant_ids)
  }
  out
}

.gsim_gbits_info <- function(haplotypes) {
  .Call(C_gsim_gbits_info, haplotypes)
}

.gsim_gbits_word <- function(haplotypes, marker, word) {
  marker <- .gsim_hapnest_integer_scalar(marker, "marker", 1)
  word <- .gsim_hapnest_integer_scalar(word, "word", 1)
  .Call(C_gsim_gbits_word, haplotypes, marker - 1L, word - 1L)
}

.gsim_gbits_copy_interval <- function(destination, destination_individual,
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
    C_gsim_gbits_copy_interval, destination, values[[1L]] - 1L,
    source, values[[2L]] - 1L, values[[3L]] - 1L, values[[4L]] - 1L
  ))
}

.gsim_gbits_copy_filtered <- function(destination, destination_individual,
                                      source, source_individual, first, last,
                                      coalescent_age, mutation_age) {
  if (!is.numeric(coalescent_age) || length(coalescent_age) != 1L ||
      !is.finite(coalescent_age) || coalescent_age < 0) {
    .gsim_stop("coalescent_age must be one finite nonnegative value.")
  }
  mutation_age <- as.double(mutation_age)
  invisible(.Call(
    C_gsim_gbits_copy_filtered, destination,
    as.integer(destination_individual - 1L), source,
    as.integer(source_individual - 1L), as.integer(first - 1L),
    as.integer(last - 1L), as.double(coalescent_age), mutation_age
  ))
}

.gsim_gbits_make_gamete <- function(destination, destination_individual,
                                    parent_h1, parent_h2, parent_individual,
                                    starting_haplotype, boundaries) {
  invisible(.Call(
    C_gsim_gbits_make_gamete, destination,
    as.integer(destination_individual - 1L), parent_h1, parent_h2,
    as.integer(parent_individual - 1L), as.integer(starting_haplotype),
    as.integer(boundaries)
  ))
}

.gsim_gbits_decode_genotypes <- function(h1, h2) {
  out <- .Call(C_gsim_gbits_decode_genotypes, h1, h2)
  dimnames(out) <- list(
    attr(h1, "sample_ids", exact = TRUE),
    attr(h1, "variant_ids", exact = TRUE)
  )
  out
}

.gsim_hapnest_founders_packed_chromosome <- function(
  backend,
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
  chromosome,
  donor_phase = "hapnest",
  return_genotypes = FALSE,
  return_segments = TRUE,
  individual_offset = 0L
) {
  chromosome <- as.character(chromosome)
  if (!length(chromosome) || anyNA(chromosome) || any(!nzchar(chromosome)) ||
      length(unique(chromosome)) != 1L) {
    .gsim_stop("The packed founder call requires exactly one chromosome label.")
  }
  plan <- .gsim_hapnest_founders(
    reference_haplotypes_h1, reference_haplotypes_h2, donor_population,
    ancestry_weights, N, Ne, rho, genetic_position, mutation_age, n, seed,
    chromosome, donor_phase, return_genotypes = FALSE,
    return_segments = TRUE, return_haplotypes = FALSE,
    individual_offset = individual_offset
  )
  reference_h1 <- .gsim_gbits_pack(backend, reference_haplotypes_h1)
  reference_h2 <- .gsim_gbits_pack(backend, reference_haplotypes_h2)
  variants <- colnames(.gsim_hapnest_raw_matrix(
    reference_haplotypes_h1, "reference_haplotypes_h1"
  ))
  ids <- paste0("syn", individual_offset + seq_len(n))
  h1 <- .gsim_gbits_zero(backend, n, length(genetic_position), ids, variants)
  h2 <- .gsim_gbits_zero(backend, n, length(genetic_position), ids, variants)
  for (i in seq_len(nrow(plan$segments))) {
    segment <- plan$segments[i, ]
    destination <- if (segment$phase == 1L) h1 else h2
    source <- if (segment$phase == 1L) reference_h1 else reference_h2
    .gsim_gbits_copy_filtered(
      destination,
      segment$individual - individual_offset,
      source,
      segment$donor_individual,
      segment$start,
      segment$end,
      segment$coalescent_age,
      mutation_age
    )
  }
  genotypes <- if (return_genotypes) .gsim_gbits_decode_genotypes(h1, h2) else NULL
  out <- list(
    h1 = h1, h2 = h2, genotypes = genotypes,
    segments = if (return_segments) plan$segments else NULL,
    sample_ids = ids, variant_ids = variants,
    chromosome = chromosome[[1L]],
    settings = c(plan$settings, list(
      storage = "gbits marker-major one-bit phased haplotypes",
      gbits_version = attr(backend, "gbits_version", exact = TRUE),
      decoded_genotypes = return_genotypes
    )),
    memory = list(
      reference_raw_bytes = 2 * length(reference_haplotypes_h1),
      reference_packed_bytes = unname(.gsim_gbits_info(reference_h1)[[4L]] +
                                        .gsim_gbits_info(reference_h2)[[4L]]),
      generated_raw_bytes = 2 * n * length(genetic_position),
      generated_packed_bytes = unname(.gsim_gbits_info(h1)[[4L]] +
                                        .gsim_gbits_info(h2)[[4L]]),
      event_record_bytes = as.numeric(object.size(plan$segments)),
      chromosome_temporary_bytes = 0,
      decoded_genotype_bytes = if (return_genotypes) length(genotypes) else 0
    )
  )
  class(out) <- c("gsim_hapnest_founders_packed", "list")
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
  founder_h1 <- .gsim_hapnest_raw_matrix(
    founder_haplotypes$h1, "founder_haplotypes$h1"
  )
  founder_h2 <- .gsim_hapnest_raw_matrix(
    founder_haplotypes$h2, "founder_haplotypes$h2"
  )
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
  founder_ids <- canonical[founder]
  supplied_ids <- rownames(founder_h1)
  missing_founders <- setdiff(founder_ids, supplied_ids)
  extra_founders <- setdiff(supplied_ids, founder_ids)
  if (length(missing_founders) || length(extra_founders) ||
      nrow(founder_h1) != length(founder_ids)) {
    detail <- c(
      if (length(missing_founders)) paste0("missing: ", paste(missing_founders, collapse = ", ")),
      if (length(extra_founders)) paste0("extra: ", paste(extra_founders, collapse = ", "))
    )
    .gsim_stop("Founder haplotype IDs must exactly match pedigree founders",
               if (length(detail)) paste0(" (", paste(detail, collapse = "; "), ")") else "",
               ".")
  }
  variant_ids <- colnames(founder_h1)
  marker_count <- ncol(founder_h1)
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
  source_h1 <- .gsim_gbits_pack(backend, founder_h1)
  source_h2 <- .gsim_gbits_pack(backend, founder_h2)
  h1 <- .gsim_gbits_zero(backend, animal_count, marker_count,
                         canonical, variant_ids)
  h2 <- .gsim_gbits_zero(backend, animal_count, marker_count,
                         canonical, variant_ids)
  founder_order <- match(founder_ids, supplied_ids)
  founder_positions <- match(founder_ids, canonical)
  for (i in seq_along(founder_ids)) {
    .gsim_gbits_copy_interval(h1, founder_positions[[i]], source_h1,
                              founder_order[[i]], 1L, marker_count)
    .gsim_gbits_copy_interval(h2, founder_positions[[i]], source_h2,
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
        .gsim_gbits_make_gamete(
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
  genotypes <- if (return_genotypes) .gsim_gbits_decode_genotypes(h1, h2) else NULL
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
      storage = "gbits marker-major one-bit phased haplotypes",
      gbits_version = attr(backend, "gbits_version", exact = TRUE),
      decoded_genotypes = return_genotypes
    ),
    memory = list(
      founder_raw_bytes = 2 * length(founder_h1),
      founder_packed_bytes = unname(.gsim_gbits_info(source_h1)[[4L]] +
                                      .gsim_gbits_info(source_h2)[[4L]]),
      generated_raw_bytes = 2 * animal_count * marker_count,
      generated_packed_bytes = unname(.gsim_gbits_info(h1)[[4L]] +
                                        .gsim_gbits_info(h2)[[4L]]),
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
