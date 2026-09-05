# Experimental transactional BED orchestration over chromosome-local packed
# phases. Persistent BIM/FAM and allele-orientation metadata remain deferred.

.gsim_bed_sink_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    .gsim_stop(name, " must be TRUE or FALSE.")
  }
  value
}

.gsim_bed_sink_create <- function(
  backend,
  path,
  sample_ids,
  overwrite = FALSE,
  buffer_variants = 64L,
  provenance = list()
) {
  if (!inherits(backend, "gsim_gbits_backend")) {
    .gsim_stop("backend must be created by .gsim_gbits_backend().")
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path)) {
    .gsim_stop("path must be one nonempty string.")
  }
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  path <- file.path(parent, basename(path))
  sample_ids <- as.character(sample_ids)
  if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    .gsim_stop("sample_ids must be nonempty, unique, nonmissing identifiers.")
  }
  overwrite <- .gsim_bed_sink_flag(overwrite, "overwrite")
  buffer_variants <- .gsim_hapnest_integer_scalar(
    buffer_variants, "buffer_variants", 1
  )
  if (buffer_variants > 4096L) {
    .gsim_stop("buffer_variants must not exceed 4096.")
  }
  if (!is.list(provenance)) {
    .gsim_stop("provenance must be a list.")
  }
  state <- new.env(parent = emptyenv())
  state$pointer <- .Call(
    C_gsim_gbits_bed_sink_create, backend, enc2utf8(path),
    length(sample_ids), overwrite, buffer_variants
  )
  state$backend <- backend
  state$path <- path
  state$sample_ids <- sample_ids
  state$chromosome <- character()
  state$marker_count <- integer()
  state$start_variant <- integer()
  state$end_variant <- integer()
  state$variant_ids <- character()
  state$provenance <- provenance
  state$finalized <- FALSE
  state$cancelled <- FALSE
  class(state) <- "gsim_bed_sink"
  state
}

.gsim_bed_sink_append <- function(
  sink,
  chromosome,
  h1,
  h2,
  variant_ids = attr(h1, "variant_ids", exact = TRUE)
) {
  if (!inherits(sink, "gsim_bed_sink") || !is.environment(sink) ||
      is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental gsim BED sink.")
  }
  if (isTRUE(sink$finalized) || isTRUE(sink$cancelled)) {
    .gsim_stop("Cannot append after BED sink finalization or cancellation.")
  }
  if (!inherits(h1, "gsim_gbits_haplotypes") ||
      !inherits(h2, "gsim_gbits_haplotypes")) {
    .gsim_stop("h1 and h2 must be packed gbits haplotype handles.")
  }
  chromosome <- as.character(chromosome)
  if (length(chromosome) != 1L || is.na(chromosome) || !nzchar(chromosome)) {
    .gsim_stop("chromosome must be one nonempty label.")
  }
  if (chromosome %in% sink$chromosome) {
    .gsim_stop("Each chromosome may be appended only once.")
  }
  info_h1 <- .gsim_gbits_info(h1)
  info_h2 <- .gsim_gbits_info(h2)
  if (!identical(unname(info_h1[1:2]), unname(info_h2[1:2]))) {
    .gsim_stop("Packed H1 and H2 dimensions must be identical.")
  }
  h1_samples <- attr(h1, "sample_ids", exact = TRUE)
  h2_samples <- attr(h2, "sample_ids", exact = TRUE)
  if (!identical(h1_samples, sink$sample_ids) ||
      !identical(h2_samples, sink$sample_ids)) {
    .gsim_stop("Packed H1/H2 sample order must exactly match sink sample_ids.")
  }
  h1_variants <- attr(h1, "variant_ids", exact = TRUE)
  h2_variants <- attr(h2, "variant_ids", exact = TRUE)
  variant_ids <- as.character(variant_ids)
  marker_count <- as.integer(info_h1[[2L]])
  if (length(variant_ids) != marker_count || anyNA(variant_ids) ||
      any(!nzchar(variant_ids)) || anyDuplicated(variant_ids) ||
      !identical(h1_variants, variant_ids) ||
      !identical(h2_variants, variant_ids)) {
    .gsim_stop("variant_ids must exactly match unique packed H1/H2 variant order.")
  }
  if (any(variant_ids %in% sink$variant_ids)) {
    .gsim_stop("Variant identifiers must be unique across appended chromosomes.")
  }

  invisible(.Call(C_gsim_gbits_bed_sink_append, sink$pointer, h1, h2))
  start <- length(sink$variant_ids) + 1L
  sink$chromosome <- c(sink$chromosome, chromosome)
  sink$marker_count <- c(sink$marker_count, marker_count)
  sink$start_variant <- c(sink$start_variant, start)
  sink$end_variant <- c(sink$end_variant, start + marker_count - 1L)
  sink$variant_ids <- c(sink$variant_ids, variant_ids)
  invisible(sink)
}

.gsim_bed_sink_info <- function(sink) {
  if (!inherits(sink, "gsim_bed_sink") || is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental gsim BED sink.")
  }
  .Call(C_gsim_gbits_bed_sink_info, sink$pointer)
}

.gsim_bed_sink_finalize <- function(sink) {
  if (!inherits(sink, "gsim_bed_sink") || is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental gsim BED sink.")
  }
  if (isTRUE(sink$finalized)) {
    .gsim_stop("BED sink has already been finalized.")
  }
  if (isTRUE(sink$cancelled)) {
    .gsim_stop("A cancelled BED sink cannot be finalized.")
  }
  invisible(.Call(C_gsim_gbits_bed_sink_finalize, sink$pointer))
  sink$finalized <- TRUE
  native <- .gsim_bed_sink_info(sink)
  expected_bytes <- 3 + length(sink$variant_ids) *
    ceiling(length(sink$sample_ids) / 4)
  manifest <- list(
    path = sink$path,
    individual_count = length(sink$sample_ids),
    variant_count = length(sink$variant_ids),
    bytes_written = unname(native[["bytes_written"]]),
    expected_bytes = expected_bytes,
    conversion_buffer_bytes = unname(native[["conversion_buffer_bytes"]]),
    lifecycle_object_bytes = unname(native[["lifecycle_object_bytes"]]),
    state = attr(native, "state", exact = TRUE),
    sample_ids = sink$sample_ids,
    variant_ids = sink$variant_ids,
    chromosome_order = sink$chromosome,
    chromosomes = data.frame(
      chromosome = sink$chromosome,
      marker_count = sink$marker_count,
      start_variant = sink$start_variant,
      end_variant = sink$end_variant,
      stringsAsFactors = FALSE
    ),
    backend = list(
      version = attr(sink$backend, "gbits_version", exact = TRUE),
      abi = attr(sink$backend, "gbits_abi", exact = TRUE)
    ),
    provenance = sink$provenance,
    metadata = list(
      bed_mode = "PLINK 1 SNP-major",
      dosage = "H1 + H2 alternative-allele count encoded as gbits raw dosage",
      bim = "deferred; allele-1 orientation must be supplied explicitly",
      fam = "deferred"
    )
  )
  class(manifest) <- c("gsim_bed_manifest", "list")
  manifest
}

.gsim_bed_sink_cancel <- function(sink) {
  if (!inherits(sink, "gsim_bed_sink") || is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental gsim BED sink.")
  }
  if (isTRUE(sink$finalized)) {
    .gsim_stop("A finalized BED sink cannot be cancelled.")
  }
  if (isTRUE(sink$cancelled)) return(invisible(sink))
  invisible(.Call(C_gsim_gbits_bed_sink_cancel, sink$pointer))
  sink$cancelled <- TRUE
  invisible(sink)
}

# Bounded qualification helper. This deliberately materializes decoded calls
# and is not used by the production sink.
.gsim_gbits_bed_read_all <- function(backend, path, individuals, variants,
                                     sample_ids = NULL, variant_ids = NULL) {
  individuals <- .gsim_hapnest_integer_scalar(individuals, "individuals", 1)
  variants <- .gsim_hapnest_integer_scalar(variants, "variants", 1)
  out <- .Call(
    C_gsim_gbits_bed_read_all, backend, enc2utf8(path), individuals, variants
  )
  dimnames(out) <- list(sample_ids, variant_ids)
  out
}
