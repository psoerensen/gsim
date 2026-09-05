.gsim_reference_descriptor <- function(reader, provenance = list()) {
  inspected <- .gsim_hap_dataset_inspect(reader)
  out <- list(
    prefix = sub("\\.hap$", "", reader$paths[["hap"]], ignore.case = TRUE),
    paths = reader$paths,
    individual_count = inspected$individual_count,
    marker_count = inspected$marker_count,
    chromosomes = inspected$chromosomes,
    sample_ids = inspected$sample_ids,
    variant_ids = inspected$variant_ids,
    format = inspected$format,
    allele_orientation = inspected$allele_orientation,
    provenance = provenance
  )
  class(out) <- c("gsim_reference", "list")
  out
}

#' Import a phased biallelic VCF reference panel
#'
#' Imports a plain, gzip, or BGZF VCF directly into chromosome-wise
#' packed HAP v1 storage with aligned BIM and FAM metadata. The supported VCF
#' subset is deliberately bounded: retained selected-sample calls are diploid,
#' complete, phased `0|0`, `0|1`, `1|0`, or `1|1`, and retained variants are
#' biallelic uppercase A/C/G/T SNPs. `GT` may occur anywhere in FORMAT.
#' Unsupported biological records are counted and either skipped or rejected;
#' malformed VCF is always rejected. Nothing is phased, imputed, normalized,
#' flipped, split, or modified.
#'
#' @param vcf Path to a plain `.vcf`, gzip VCF, or BGZF VCF. Compression is
#'   detected from file bytes rather than the name. BCF is unsupported.
#' @param map A data frame containing `chromosome`, `genetic_position_cm`, and
#'   exactly one alignment key. `variant_id` requests exact per-variant
#'   alignment. `base_pair_position` supplies at least two strictly increasing
#'   sparse knots per imported chromosome; cumulative cM is copied at knots and
#'   linearly interpolated between them without extrapolation.
#' @param output Extension-free output prefix for `.hap`, `.bim`, and `.fam`.
#' @param sample_metadata Optional data frame keyed by exact `individual_id`.
#'   It may supply `family_id` and sex codes 0/1/2. Parent IDs must be 0. The
#'   defaults are FID `reference`, unknown sex, and missing phenotype.
#' @param samples `NULL` for all VCF samples, or unique exact sample IDs in the
#'   desired HAP/FAM order. Unselected sample genotype fields are not decoded.
#' @param chromosome `NULL` or one exact chromosome label.
#' @param region `NULL` or inclusive positive base-pair bounds `c(start, end)`.
#'   A region requires `chromosome`; the compressed stream is scanned rather
#'   than accessed through a tabix index.
#' @param unsupported Either `"skip"` to count and omit ordinary unsupported
#'   biological records or `"error"` to stop at the first one.
#' @param overwrite Whether to replace an existing complete triplet.
#'
#' @return A lightweight `gsim_reference` descriptor. Biological alleles remain
#'   in packed files and are not returned as matrices.
#' @importFrom utils head object.size
#' @export
#'
#' @details VCF REF is stored as bit 0 and BIM A2; VCF ALT is bit 1 and BIM A1.
#' The left and right phased GT alleles become H1 and H2. Memory is bounded by
#' one chromosome's two one-bit planes, one logical VCF line, a fixed 64 KiB
#' decompression buffer, and a two-byte-per-selected-sample allele buffer.
#' Sparse interpolation evaluates
#' `cm_left + (cm_right - cm_left) * (bp - bp_left) / (bp_right - bp_left)`
#' in that order using R double arithmetic. Exact knot positions copy the
#' supplied value. Duplicate knots and repeated chromosome blocks are errors.
gsim_import_vcf <- function(
  vcf, map, output, sample_metadata = NULL, samples = NULL,
  chromosome = NULL, region = NULL, unsupported = c("skip", "error"),
  overwrite = FALSE
) {
  unsupported <- match.arg(unsupported)
  manifest <- .gsim_import_vcf_internal(
    vcf, map, output,
    sample_metadata = sample_metadata, samples = samples,
    chromosome = chromosome, region = region, unsupported = unsupported,
    overwrite = overwrite
  )
  reader <- .gsim_hap_dataset_open(output)
  on.exit(.gsim_hap_dataset_close(reader), add = TRUE)
  descriptor <- .gsim_reference_descriptor(reader, manifest$provenance)
  descriptor$import <- manifest$import
  descriptor$manifest <- manifest
  descriptor
}

#' Open a packed HAP/BIM/FAM reference panel
#'
#' Reconstructs a lightweight validated descriptor for an existing phased
#' dataset. HAP counts and chromosome ranges are checked against BIM/FAM.
#'
#' @param prefix Extension-free HAP/BIM/FAM prefix.
#' @return A `gsim_reference` descriptor containing paths and stable identities.
#' @export
gsim_reference <- function(prefix) {
  reader <- .gsim_hap_dataset_open(prefix)
  on.exit(.gsim_hap_dataset_close(reader), add = TRUE)
  .gsim_reference_descriptor(
    reader, list(operation = "validated existing HAP/BIM/FAM reference")
  )
}

.gsim_public_align <- function(value, ids, name, numeric = FALSE) {
  .gsim_hapnest_align_named_reference_vector(value, ids, name, numeric)
}

.gsim_public_phased_input <- function(value, name) {
  if (is.character(value) && length(value) == 1L) {
    value <- gsim_reference(value)
  }
  if (!inherits(value, "gsim_reference") || !is.character(value$prefix) ||
      length(value$prefix) != 1L || is.na(value$prefix) ||
      !nzchar(value$prefix)) {
    .gsim_stop(name, " must be a phased-population descriptor or one prefix.")
  }
  value
}

.gsim_public_threads <- function(threads) {
  .gsim_hapnest_integer_scalar(threads, "threads", 1)
}

.gsim_public_founder_ids <- function(n, founder_ids) {
  if (is.null(n) == is.null(founder_ids)) {
    .gsim_stop("Exactly one of n or founder_ids must define the base population.")
  }
  if (!is.null(n)) {
    n <- .gsim_hapnest_integer_scalar(n, "n", 1)
    return(paste0("syn", seq_len(n)))
  }
  if (!is.character(founder_ids)) {
    .gsim_stop("founder_ids must be a character vector.")
  }
  founder_ids <- enc2utf8(founder_ids)
  if (!length(founder_ids) || anyNA(founder_ids) ||
      any(!nzchar(founder_ids)) || anyDuplicated(founder_ids)) {
    .gsim_stop("founder_ids must be unique, nonmissing, and nonempty.")
  }
  founder_ids
}

.gsim_simulate_founders_impl <- function(
  reference, n = NULL, founder_ids = NULL, populations, ancestry_weights,
  mutation_age, N, Ne, rho, seed, output, batch_size = NULL, threads = 1L,
  overwrite = FALSE, return_segments = FALSE
) {
  reference <- .gsim_public_phased_input(reference, "reference")
  founder_ids <- .gsim_public_founder_ids(n, founder_ids)
  threads <- .gsim_public_threads(threads)
  requested_batch_size <- if (is.null(batch_size)) 8192L else
    .gsim_hapnest_integer_scalar(batch_size, "batch_size", 1)
  bounded_batch_size <- min(requested_batch_size, length(founder_ids))
  actual_batch_size <- if (bounded_batch_size > .Machine$integer.max - 63L) {
    length(founder_ids)
  } else {
    min(length(founder_ids),
        as.integer(((bounded_batch_size + 63L) %/% 64L) * 64L))
  }
  reader <- .gsim_hap_dataset_open(reference$prefix)
  on.exit(try(.gsim_hap_dataset_close(reader), silent = TRUE), add = TRUE)
  populations <- .gsim_public_align(
    populations, reader$samples$individual_id, "populations"
  )
  mutation_age <- .gsim_public_align(
    mutation_age, reader$variants$variant_id, "mutation_age", numeric = TRUE
  )
  names(mutation_age) <- reader$variants$variant_id
  sample_metadata <- .gsim_plink_sample_metadata(
    founder_ids, family_id = rep.int("base", length(founder_ids))
  )
  provenance <- list(
    operation = "HAPNEST-compatible synthetic base population",
    reference = reader$paths,
    founder_model = "phase-specific donor copying",
    seed = seed, founder_ids = "explicit final order"
  )
  sink <- .gsim_hap_dataset_create(
    output, sample_metadata,
    overwrite, provenance
  )
  completed <- FALSE
  on.exit({
    if (!completed) try(.gsim_hap_dataset_cancel(sink), silent = TRUE)
  }, add = TRUE)
  peak_payload <- 0
  audit <- list()
  for (chromosome in reader$chromosome) {
    index <- match(chromosome, reader$chromosome)
    start <- sum(head(reader$marker_count, index - 1L)) + 1L
    rows <- seq.int(start, length.out = reader$marker_count[[index]])
    variants <- reader$variants[rows, , drop = FALSE]
    cm <- stats::setNames(variants$genetic_position_cm, variants$variant_id)
    reference_handles <- .gsim_hap_dataset_load_chromosome(reader, chromosome)
    on.exit({
      try(.gsim_packed_close(reference_handles$h1), silent = TRUE)
      try(.gsim_packed_close(reference_handles$h2), silent = TRUE)
    }, add = TRUE)
    .gsim_hap_dataset_begin_chromosome(sink, chromosome, variants)
    chromosome_audit <- list()
    batch_number <- 0L
    starts <- seq.int(0L, length(founder_ids) - 1L, by = actual_batch_size)
    for (individual_offset in starts) {
      batch_number <- batch_number + 1L
      count <- min(actual_batch_size, length(founder_ids) - individual_offset)
      batch_ids <- founder_ids[individual_offset + seq_len(count)]
      batch_result <- local({
        batch <- .gsim_hapnest_founders_packed_reference_chromosome(
          reference_handles$h1, reference_handles$h2,
          stats::setNames(populations, reader$samples$individual_id),
          ancestry_weights, N, Ne, rho, unname(cm),
          unname(mutation_age[variants$variant_id]), count, seed, chromosome,
          return_genotypes = FALSE, return_segments = return_segments,
          individual_offset = individual_offset, threads = threads
        )
        on.exit({
          try(.gsim_packed_close(batch$h1), silent = TRUE)
          try(.gsim_packed_close(batch$h2), silent = TRUE)
        }, add = TRUE)
        batch$h1 <- .gsim_packed_tag(batch$h1, batch_ids, variants$variant_id)
        batch$h2 <- .gsim_packed_tag(batch$h2, batch_ids, variants$variant_id)
        .gsim_hap_dataset_write_batch(
          sink, batch$h1, batch$h2, individual_offset
        )
        list(segments = batch$segments, memory = batch$memory)
      })
      if (return_segments) {
        chromosome_audit[[batch_number]] <- batch_result$segments
      }
      peak_payload <- max(
        peak_payload,
        batch_result$memory$reference_packed_bytes +
          batch_result$memory$generated_packed_bytes
      )
    }
    if (return_segments) {
      audit[[chromosome]] <- do.call(rbind, chromosome_audit)
      rownames(audit[[chromosome]]) <- NULL
    }
    .gsim_packed_close(reference_handles$h1)
    .gsim_packed_close(reference_handles$h2)
    reference_handles <- list(h1 = NULL, h2 = NULL)
  }
  manifest <- .gsim_hap_dataset_finalize(sink)
  completed <- TRUE
  output_reader <- .gsim_hap_dataset_open(output)
  on.exit(.gsim_hap_dataset_close(output_reader), add = TRUE)
  result <- .gsim_reference_descriptor(output_reader, provenance)
  result$manifest <- manifest
  result$simulation <- list(
    population = "synthetic base/founder population",
    reference_prefix = reference$prefix,
    seed = seed, founder_count = length(founder_ids),
    founder_ids = founder_ids,
    requested_batch_size = if (is.null(batch_size)) NULL else batch_size,
    actual_batch_size = actual_batch_size,
    batch_alignment = "rounded up to a 64-sample packed-word boundary",
    threads = threads,
    peak_chromosome_biological_payload_bytes = peak_payload,
    dense_matrix_allocated = FALSE,
    rng = "chromosome-label-keyed SplitMix64; R global RNG unused"
  )
  if (return_segments) result$segment_audit <- audit
  class(result) <- c("gsim_base_population", class(result))
  result
}

#' Simulate a reusable phased synthetic base population
#'
#' Generates unrelated founders from a packed reference panel using the
#' HAPNEST-compatible historical copying model. Output is always phased
#' HAP/BIM/FAM so that the base population can be reused by independent
#' pedigree simulations.
#'
#' @param reference A [gsim_reference()] descriptor or HAP/BIM/FAM prefix.
#' @param n Number of founders. Exactly one of `n` and `founder_ids` is used.
#'   Generated IDs are `syn1`, `syn2`, and so on.
#' @param founder_ids Explicit unique character founder IDs in final sample
#'   order.
#' @param populations Population labels named by reference FAM IID.
#' @param ancestry_weights Named donor-population probabilities.
#' @param mutation_age Finite mutation ages named by reference BIM ID.
#' @param N,Ne,rho Named HAPNEST population parameters in the committed
#'   cumulative-cM founder convention.
#' @param seed Founder-generation seed.
#' @param output Extension-free HAP/BIM/FAM output prefix.
#' @param batch_size Positive requested founder batch size. Internal batches are
#'   aligned to packed 64-sample words; `NULL` uses a bounded default.
#' @param threads Positive native founder worker count.
#' @param overwrite Whether to replace an existing complete triplet.
#'
#' @return A lightweight reusable phased base-population descriptor.
#' @export
#'
#' @details Founder H1 copies only donor H1 and founder H2 only donor H2 under
#' the chromosome-label-keyed SplitMix64 contract and strict
#' `T < mutation_age` filtering. `batch_size = NULL` requests 8,192 founders;
#' explicit sizes are rounded up to a 64-sample word boundary and the last batch
#' may be partial. Each batch is materialized by static native workers that own
#' disjoint packed words, then written directly to its final HAP sample range.
#' Batch size and thread count do not enter RNG streams. The R global RNG is not
#' used.
gsim_simulate_founders <- function(
  reference, n = NULL, founder_ids = NULL, populations, ancestry_weights,
  mutation_age, N, Ne, rho, seed, output, batch_size = NULL, threads = 1L,
  overwrite = FALSE
) {
  .gsim_simulate_founders_impl(
    reference, n, founder_ids, populations, ancestry_weights, mutation_age,
    N, Ne, rho, seed, output, batch_size, threads, overwrite, FALSE
  )
}

.gsim_pedigree_founder_ids <- function(pedigree) {
  canonical <- as.character(pedigree$canonical_order)
  tab <- pedigree$pedigree[
    match(canonical, as.character(pedigree$pedigree$animal)), , drop = FALSE
  ]
  canonical[is.na(tab$sire) & is.na(tab$dam)]
}

.gsim_simulate_pedigree_impl <- function(
  founders, pedigree, seed, output, format = c("hap", "bed"),
  overwrite = FALSE, return_crossovers = FALSE
) {
  format <- match.arg(format)
  founders <- .gsim_public_phased_input(founders, "founders")
  if (!inherits(pedigree, "gsim_pedigree")) {
    .gsim_stop("pedigree must be a gsim_pedigree object.")
  }
  reader <- .gsim_hap_dataset_open(founders$prefix)
  on.exit(try(.gsim_hap_dataset_close(reader), silent = TRUE), add = TRUE)
  canonical <- as.character(pedigree$canonical_order)
  founder_ids <- .gsim_pedigree_founder_ids(pedigree)
  if (!length(founder_ids)) .gsim_stop("pedigree contains no founders.")
  missing <- setdiff(founder_ids, reader$samples$individual_id)
  extra <- setdiff(reader$samples$individual_id, founder_ids)
  if (length(missing) || length(extra) ||
      length(founder_ids) != nrow(reader$samples)) {
    detail <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(extra)) paste0("extra: ", paste(extra, collapse = ", "))
    )
    .gsim_stop(
      "Pedigree founders must exactly match phased base-population IDs",
      if (length(detail)) paste0(" (", paste(detail, collapse = "; "), ")") else "",
      "."
    )
  }
  sample_metadata <- .gsim_plink_pedigree_metadata(pedigree)
  provenance <- list(
    operation = "Mendelian pedigree simulation from phased base population",
    founders = reader$paths,
    pedigree_model = "chromosome-wise Poisson no-interference meiosis",
    meiosis_map_conversion = "BIM cumulative cM divided by 100 to Morgans",
    seed = seed, output_format = format
  )
  sink <- if (format == "hap") {
    .gsim_hap_dataset_create(output, sample_metadata, overwrite, provenance)
  } else {
    .gsim_plink_dataset_create(output, sample_metadata, overwrite,
                               provenance = provenance)
  }
  completed <- FALSE
  on.exit({
    if (!completed) {
      if (format == "hap") try(.gsim_hap_dataset_cancel(sink), silent = TRUE)
      else try(.gsim_plink_dataset_cancel(sink), silent = TRUE)
    }
  }, add = TRUE)
  peak_payload <- 0
  process <- function(chromosome) {
    index <- match(chromosome, reader$chromosome)
    start <- sum(head(reader$marker_count, index - 1L)) + 1L
    rows <- seq.int(start, length.out = reader$marker_count[[index]])
    variants <- reader$variants[rows, , drop = FALSE]
    base <- .gsim_hap_dataset_load_chromosome(reader, chromosome)
    on.exit({
      try(.gsim_packed_close(base$h1), silent = TRUE)
      try(.gsim_packed_close(base$h2), silent = TRUE)
    }, add = TRUE)
    descendants <- .gsim_pedigree_genotypes_packed_chromosome(
      pedigree, list(h1 = base$h1, h2 = base$h2),
      rep.int(chromosome, nrow(variants)),
      variants$genetic_position_cm / 100, seed,
      return_haplotypes = TRUE, return_genotypes = FALSE,
      return_crossovers = return_crossovers
    )
    on.exit({
      try(.gsim_packed_close(descendants$h1), silent = TRUE)
      try(.gsim_packed_close(descendants$h2), silent = TRUE)
    }, add = TRUE)
    if (format == "hap") {
      .gsim_hap_dataset_append(sink, chromosome, descendants$h1,
                               descendants$h2, variants)
    } else {
      .gsim_plink_dataset_append(sink, chromosome, descendants$h1,
                                 descendants$h2, variants)
    }
    if (return_crossovers) audit[[chromosome]] <<- descendants$crossover_audit
    descendants$memory$founder_packed_bytes +
      descendants$memory$generated_packed_bytes
  }
  audit <- list()
  for (chromosome in reader$chromosome) {
    peak_payload <- max(peak_payload, process(chromosome))
  }
  manifest <- if (format == "hap") .gsim_hap_dataset_finalize(sink) else
    .gsim_plink_dataset_finalize(sink)
  completed <- TRUE
  manifest$simulation <- list(
    population = "pedigreed population",
    founder_prefix = founders$prefix,
    founder_count = length(founder_ids), pedigree_count = length(canonical),
    peak_chromosome_biological_payload_bytes = peak_payload,
    dense_genotype_matrix_allocated = FALSE,
    genetic_position = list(stored_map = "cumulative cM",
                            meiosis = "cumulative Morgans = BIM cM / 100"),
    rng = "SplitMix64; R global RNG unused", seed = seed, threads = 1L,
    parallelism = paste(
      "single-threaded: generations are dependent and marker-major packed",
      "words can contain both parents and children"
    )
  )
  if (return_crossovers) manifest$crossover_audit <- audit
  class(manifest) <- c("gsim_simulation_manifest", class(manifest))
  manifest
}

#' Simulate Mendelian descendants from a phased base population
#'
#' Extends an existing phased founder population through chromosome-wise,
#' no-interference meiosis. Founder-model parameters are deliberately absent:
#' historical founder copying and future biological meiosis are separate models.
#'
#' @param founders A phased base-population descriptor or HAP/BIM/FAM prefix.
#'   Its FAM IDs must exactly equal the pedigree founder IDs.
#' @param pedigree A [gsim_pedigree()] object with parent-before-offspring order.
#' @param seed Pedigree-meiosis seed; it cannot alter the stored founders.
#' @param output Extension-free output prefix.
#' @param format `"hap"` for reusable phased HAP/BIM/FAM or `"bed"` for
#'   unphased dosage BED/BIM/FAM.
#' @param overwrite Whether to replace an existing complete triplet.
#'
#' @return A compact output manifest; no dense biological matrix is returned.
#' @export
#'
#' @details Stored founder H1/H2 labels are preserved. In descendants H1 is the
#' paternal gamete and H2 is the maternal gamete. BIM cumulative cM is divided
#' by 100 for biological meiosis in Morgans. Bit 1 remains ALT/BIM A1 and BED
#' dosage remains H1 + H2. Pedigree execution is currently single-threaded
#' because dependency-safe generations can still share mutable packed words
#' with their parents.
gsim_simulate_pedigree <- function(
  founders, pedigree, seed, output, format = c("hap", "bed"),
  overwrite = FALSE
) {
  .gsim_simulate_pedigree_impl(
    founders, pedigree, seed, output, format, overwrite, FALSE
  )
}
