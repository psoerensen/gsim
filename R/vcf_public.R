.gsim_public_backends <- function() {
  list(packed = .gsim_packed_backend(), metadata = .gsim_metadata_backend())
}

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
#' Imports an ordinary uncompressed text VCF directly into chromosome-wise
#' packed HAP v1 storage with aligned BIM and FAM metadata. The supported VCF
#' subset is deliberately strict: samples must be diploid, every `GT` must be
#' phased and one of `0|0`, `0|1`, `1|0`, or `1|1`, and variants must be
#' biallelic uppercase A/C/G/T SNPs without missing calls. `GT` may occur at
#' any position in FORMAT. No records are phased, imputed, normalized, flipped,
#' split, or discarded.
#'
#' @param vcf Path to an uncompressed `.vcf` file.
#' @param map A data frame containing `chromosome`, `genetic_position_cm`, and
#'   exactly one alignment key: `variant_id` or `base_pair_position`. Every VCF
#'   variant must match exactly once. Positions are cumulative cM and are not
#'   interpolated or converted during import.
#' @param output Extension-free output prefix for `.hap`, `.bim`, and `.fam`.
#' @param sample_metadata Optional data frame keyed by exact `individual_id`.
#'   It may supply `family_id` and sex codes 0/1/2. Parent IDs must be 0. The
#'   defaults are FID `reference`, unknown sex, and missing phenotype.
#' @param overwrite Whether to replace an existing complete triplet.
#'
#' @return A lightweight `gsim_reference` descriptor. Biological alleles remain
#'   in packed files and are not returned as matrices.
#' @importFrom utils head object.size
#' @export
#'
#' @details VCF REF is stored as bit 0 and BIM A2; VCF ALT is bit 1 and BIM A1.
#' The left and right phased GT alleles become H1 and H2. Memory is bounded by
#' one chromosome's two one-bit planes plus a two-byte-per-sample record buffer.
gsim_import_vcf <- function(vcf, map, output, sample_metadata = NULL,
                            overwrite = FALSE) {
  backend <- .gsim_public_backends()
  manifest <- .gsim_import_vcf_internal(
    backend$packed, backend$metadata, vcf, map, output, sample_metadata, overwrite
  )
  reader <- .gsim_hap_dataset_open(backend$packed, backend$metadata, output)
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
  backend <- .gsim_public_backends()
  reader <- .gsim_hap_dataset_open(backend$packed, backend$metadata, prefix)
  on.exit(.gsim_hap_dataset_close(reader), add = TRUE)
  .gsim_reference_descriptor(
    reader, list(operation = "validated existing HAP/BIM/FAM reference")
  )
}

.gsim_public_align <- function(value, ids, name, numeric = FALSE) {
  .gsim_hapnest_align_named_reference_vector(value, ids, name, numeric)
}

#' Simulate packed founders and Mendelian pedigree descendants
#'
#' Runs the qualified HAPNEST-compatible founder generator and marker-level
#' no-interference pedigree meiosis directly from a packed reference dataset.
#' Each chromosome is loaded, simulated, written, and released before the next.
#'
#' @param reference A `gsim_reference` descriptor or existing HAP/BIM/FAM prefix.
#' @param pedigree A [gsim_pedigree()] object. Every individual with both parents
#'   missing is a founder; one-known-parent records remain unsupported.
#' @param populations Population labels named by every reference FAM IID.
#' @param ancestry_weights Named donor-population probabilities.
#' @param mutation_age Finite mutation ages named by every reference BIM ID.
#' @param N,Ne,rho Named HAPNEST population parameters. Genetic positions and
#'   `rho` retain the committed founder-core cumulative-cM convention.
#' @param seed Nonnegative exactly represented integer seed.
#' @param output Extension-free output prefix.
#' @param format Either `"hap"` for phased HAP/BIM/FAM or `"bed"` for
#'   dosage BED/BIM/FAM.
#' @param overwrite Whether to replace an existing complete output triplet.
#'
#' @return A compact simulation manifest; no dense genotype matrix is returned.
#' @export
#'
#' @details Founder H1 copies only donor H1 and founder H2 only donor H2, using
#' the unchanged chromosome-label-keyed SplitMix64 event plan and strict
#' `T < mutation_age` filter. In descendants H1 is the paternal gamete and H2
#' the maternal gamete. Biological meiosis requires Morgans, so this public
#' workflow explicitly converts BIM cumulative cM to Morgans by division by
#' 100. Bit 1 remains ALT/BIM A1 and BED dosage remains H1 + H2. The R global
#' RNG is never used.
gsim_simulate <- function(
  reference, pedigree, populations, ancestry_weights, mutation_age,
  N, Ne, rho, seed, output, format = c("hap", "bed"), overwrite = FALSE
) {
  format <- match.arg(format)
  if (is.character(reference) && length(reference) == 1L) {
    reference <- gsim_reference(reference)
  }
  if (!inherits(reference, "gsim_reference") ||
      !is.character(reference$prefix) || length(reference$prefix) != 1L) {
    .gsim_stop("reference must be a gsim_reference descriptor or one prefix.")
  }
  if (!inherits(pedigree, "gsim_pedigree")) {
    .gsim_stop("pedigree must be a gsim_pedigree object.")
  }
  backend <- .gsim_public_backends()
  reader <- .gsim_hap_dataset_open(
    backend$packed, backend$metadata, reference$prefix
  )
  on.exit(try(.gsim_hap_dataset_close(reader), silent = TRUE), add = TRUE)
  populations <- .gsim_public_align(
    populations, reader$samples$individual_id, "populations"
  )
  mutation_age <- .gsim_public_align(
    mutation_age, reader$variants$variant_id, "mutation_age", numeric = TRUE
  )
  names(mutation_age) <- reader$variants$variant_id
  canonical <- as.character(pedigree$canonical_order)
  tab <- pedigree$pedigree[
    match(canonical, as.character(pedigree$pedigree$animal)), , drop = FALSE]
  founder_ids <- canonical[is.na(tab$sire) & is.na(tab$dam)]
  if (!length(founder_ids)) .gsim_stop("pedigree contains no founders.")
  sample_metadata <- .gsim_plink_pedigree_metadata(pedigree)
  provenance <- list(
    operation = "packed founder and pedigree simulation",
    reference = reader$paths,
    founder_model = "HAPNEST-compatible phase-specific donor copying",
    pedigree_model = "chromosome-wise Poisson no-interference meiosis",
    meiosis_map_conversion = "BIM cumulative cM divided by 100 to Morgans",
    seed = seed, output_format = format
  )
  sink <- if (format == "hap") {
    .gsim_hap_dataset_create(backend$packed, backend$metadata, output,
                             sample_metadata, overwrite, provenance)
  } else {
    .gsim_plink_dataset_create(backend$packed, backend$metadata, output,
                               sample_metadata, overwrite,
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
    cm <- stats::setNames(variants$genetic_position_cm, variants$variant_id)
    ages <- mutation_age[variants$variant_id]
    founders <- .gsim_hapnest_founders_from_hap_chromosome(
      reader, chromosome, stats::setNames(populations,
                                           reader$samples$individual_id),
      ancestry_weights, N, Ne, rho, cm, ages, length(founder_ids), seed,
      return_genotypes = FALSE, return_segments = FALSE
    )
    on.exit({
      try(.gsim_packed_close(founders$h1), silent = TRUE)
      try(.gsim_packed_close(founders$h2), silent = TRUE)
    }, add = TRUE)
    founders$h1 <- .gsim_packed_tag(founders$h1, founder_ids,
                                    variants$variant_id)
    founders$h2 <- .gsim_packed_tag(founders$h2, founder_ids,
                                    variants$variant_id)
    descendants <- .gsim_pedigree_genotypes_packed_chromosome(
      backend$packed, pedigree, list(h1 = founders$h1, h2 = founders$h2),
      rep.int(chromosome, nrow(variants)),
      variants$genetic_position_cm / 100, seed,
      return_haplotypes = TRUE, return_genotypes = FALSE,
      return_crossovers = FALSE
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
    founders$memory$reference_packed_bytes +
      founders$memory$generated_packed_bytes +
      descendants$memory$generated_packed_bytes
  }
  for (chromosome in reader$chromosome) {
    peak_payload <- max(peak_payload, process(chromosome))
  }
  manifest <- if (format == "hap") .gsim_hap_dataset_finalize(sink) else
    .gsim_plink_dataset_finalize(sink)
  completed <- TRUE
  manifest$simulation <- list(
    reference_prefix = reference$prefix,
    founder_count = length(founder_ids), pedigree_count = length(canonical),
    peak_chromosome_biological_payload_bytes = peak_payload,
    dense_genotype_matrix_allocated = FALSE,
    genetic_position = list(founder = "cumulative cM",
                            meiosis = "cumulative Morgans = BIM cM / 100"),
    rng = "SplitMix64; R global RNG unused"
  )
  class(manifest) <- c("gsim_simulation_manifest", class(manifest))
  manifest
}
