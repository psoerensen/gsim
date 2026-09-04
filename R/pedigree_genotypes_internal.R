# Internal byte-matrix oracle for marker-level Mendelian inheritance.  The
# frozen model and ownership boundary are in
# docs/design/pedigree_marker_meiosis.md.  This API is intentionally unexported.

.gsim_meiosis_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .gsim_stop(name, " must be TRUE or FALSE.")
  }
  x
}

.gsim_meiosis_raw_vector <- function(x, name) {
  if (is.matrix(x) || !length(x) ||
      !(is.raw(x) || is.logical(x) || is.integer(x) || is.double(x))) {
    .gsim_stop(name, " must be a nonempty 0/1 vector.")
  }
  if (!is.raw(x) && (anyNA(x) || any(!is.finite(x)) ||
                     any(!(x %in% c(0, 1))))) {
    .gsim_stop(name, " must contain only nonmissing 0/1 alleles.")
  }
  if (is.raw(x) && any(!(as.integer(x) %in% c(0L, 1L)))) {
    .gsim_stop(name, " must contain only 0/1 alleles.")
  }
  as.raw(x)
}

.gsim_meiosis_materialize <- function(
  parent_h1,
  parent_h2,
  genetic_position,
  crossovers = numeric(),
  starting_haplotype = 1L
) {
  parent_h1 <- .gsim_meiosis_raw_vector(parent_h1, "parent_h1")
  parent_h2 <- .gsim_meiosis_raw_vector(parent_h2, "parent_h2")
  if (length(parent_h1) != length(parent_h2)) {
    .gsim_stop("parent_h1 and parent_h2 must have identical lengths.")
  }
  genetic_position <- as.double(genetic_position)
  if (length(genetic_position) != length(parent_h1) ||
      any(!is.finite(genetic_position)) ||
      any(diff(genetic_position) < 0)) {
    .gsim_stop(
      "genetic_position must be finite, nondecreasing, and aligned to the parental haplotypes."
    )
  }
  crossovers <- as.double(crossovers)
  if (any(!is.finite(crossovers)) || any(diff(crossovers) < 0) ||
      any(crossovers < genetic_position[[1L]]) ||
      any(crossovers > genetic_position[[length(genetic_position)]])) {
    .gsim_stop(
      "crossovers must be finite, nondecreasing, and inside the chromosome interval."
    )
  }
  starting_haplotype <- .gsim_hapnest_integer_scalar(
    starting_haplotype, "starting_haplotype", 1
  )
  if (starting_haplotype > 2L) {
    .gsim_stop("starting_haplotype must be 1 or 2.")
  }
  .Call(
    C_gsim_meiosis_materialize,
    parent_h1,
    parent_h2,
    genetic_position,
    crossovers,
    starting_haplotype
  )
}

.gsim_meiosis_empty_meioses <- function() {
  data.frame(
    animal = character(), parent = character(), parental_side = character(),
    child_phase = character(), chromosome = character(),
    starting_haplotype = integer(), crossover_count = integer(),
    stringsAsFactors = FALSE
  )
}

.gsim_meiosis_empty_crossovers <- function() {
  data.frame(
    animal = character(), parent = character(), parental_side = character(),
    child_phase = character(), chromosome = character(),
    crossover_index = integer(), genetic_position = numeric(),
    stringsAsFactors = FALSE
  )
}

.gsim_pedigree_genotypes <- function(
  pedigree,
  founder_haplotypes,
  chromosome,
  genetic_position,
  seed,
  return_haplotypes = TRUE,
  return_genotypes = TRUE,
  return_crossovers = TRUE,
  batch_size = NULL
) {
  if (!inherits(pedigree, "gsim_pedigree") ||
      !is.data.frame(pedigree$pedigree) ||
      is.null(pedigree$canonical_order)) {
    .gsim_stop("pedigree must be a gsim_pedigree object.")
  }
  canonical <- as.character(pedigree$canonical_order)
  tab <- pedigree$pedigree
  required_columns <- c("animal", "sire", "dam")
  if (any(!required_columns %in% names(tab)) || !length(canonical) ||
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
    .gsim_stop(
      "Animal '", tab$animal[[i]], "' has one known parent; missing ",
      missing_side, " is not supported."
    )
  }
  for (i in which(!founder)) {
    if (is.na(sire_position[[i]]) || is.na(dam_position[[i]])) {
      .gsim_stop("Recorded parent of animal '", tab$animal[[i]],
                 "' is absent from canonical_order.")
    }
    if (sire_position[[i]] >= i || dam_position[[i]] >= i) {
      .gsim_stop(
        "Canonical parent-before-offspring contract is violated for animal '",
        tab$animal[[i]], "'."
      )
    }
  }

  if (!is.list(founder_haplotypes) || is.null(names(founder_haplotypes)) ||
      anyDuplicated(names(founder_haplotypes)) ||
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
  genetic_position <- as.double(genetic_position)
  if (length(chromosome) != marker_count || anyNA(chromosome) ||
      any(!nzchar(chromosome))) {
    .gsim_stop("chromosome must provide one nonmissing label per variant.")
  }
  if (length(genetic_position) != marker_count ||
      any(!is.finite(genetic_position))) {
    .gsim_stop("genetic_position must provide one finite value per variant.")
  }
  runs <- rle(chromosome)
  if (anyDuplicated(runs$values)) {
    .gsim_stop("Each chromosome must occupy one contiguous variant block.")
  }
  block_end <- cumsum(runs$lengths)
  block_start <- c(1L, head(block_end, -1L) + 1L)
  for (block in seq_along(block_start)) {
    at <- block_start[[block]]:block_end[[block]]
    if (length(at) > 1L && any(diff(genetic_position[at]) < 0)) {
      .gsim_stop("genetic_position must be nondecreasing within each chromosome.")
    }
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
  h1 <- matrix(as.raw(0), nrow = animal_count, ncol = marker_count,
               dimnames = list(canonical, variant_ids))
  h2 <- h1
  founder_order <- match(founder_ids, supplied_ids)
  h1[match(founder_ids, canonical), ] <- founder_h1[founder_order, , drop = FALSE]
  h2[match(founder_ids, canonical), ] <- founder_h2[founder_order, , drop = FALSE]

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
        side_name <- if (side == 1L) "paternal" else "maternal"
        child_phase <- if (side == 1L) "H1" else "H2"
        for (block in seq_along(block_start)) {
          at <- block_start[[block]]:block_end[[block]]
          draw <- .Call(
            C_gsim_meiosis_draw,
            h1[parent_position, at],
            h2[parent_position, at],
            genetic_position[at],
            as.double(seed),
            canonical[[i]],
            runs$values[[block]],
            as.integer(side),
            return_crossovers
          )
          gamete <- if (return_crossovers) draw$gamete else draw
          if (side == 1L) {
            h1[i, at] <- gamete
          } else {
            h2[i, at] <- gamete
          }
          if (return_crossovers) {
            meiosis_number <- meiosis_number + 1L
            meiosis_rows[[meiosis_number]] <- data.frame(
              animal = canonical[[i]], parent = parent_id,
              parental_side = side_name, child_phase = child_phase,
              chromosome = runs$values[[block]],
              starting_haplotype = draw$starting_haplotype,
              crossover_count = length(draw$crossovers),
              stringsAsFactors = FALSE
            )
            if (length(draw$crossovers)) {
              crossover_number <- crossover_number + 1L
              crossover_rows[[crossover_number]] <- data.frame(
                animal = rep.int(canonical[[i]], length(draw$crossovers)),
                parent = rep.int(parent_id, length(draw$crossovers)),
                parental_side = rep.int(side_name, length(draw$crossovers)),
                child_phase = rep.int(child_phase, length(draw$crossovers)),
                chromosome = rep.int(runs$values[[block]], length(draw$crossovers)),
                crossover_index = seq_along(draw$crossovers),
                genetic_position = draw$crossovers,
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
  }

  genotype <- if (return_genotypes) .gsim_hapnest_pair(h1, h2) else NULL
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
  chromosome_blocks <- data.frame(
    chromosome = runs$values,
    start_variant = block_start,
    end_variant = block_end,
    length_morgans = genetic_position[block_end] - genetic_position[block_start],
    stringsAsFactors = FALSE
  )
  pedigree_alignment <- data.frame(
    animal = canonical,
    canonical_index = seq_along(canonical),
    external_index = match(canonical, as.character(pedigree$pedigree$animal)),
    sire = tab$sire,
    dam = tab$dam,
    founder = founder,
    stringsAsFactors = FALSE
  )
  variant_map <- data.frame(
    variant = variant_ids,
    variant_index = seq_along(variant_ids),
    chromosome = chromosome,
    genetic_position = genetic_position,
    stringsAsFactors = FALSE
  )
  out <- list(
    h1 = if (return_haplotypes) h1 else NULL,
    h2 = if (return_haplotypes) h2 else NULL,
    genotypes = genotype,
    crossover_audit = crossover_audit,
    sample_ids = canonical,
    variant_ids = variant_ids,
    variant_map = variant_map,
    chromosome_blocks = chromosome_blocks,
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
      seed = as.double(seed),
      global_rng = "not used",
      batch_policy = "operational batch boundaries do not enter RNG streams",
      unrelated_animal_policy = "stable IDs leave existing streams unchanged"
    )
  )
  class(out) <- c("gsim_pedigree_genotypes", "list")
  out
}
