# Experimental chromosome-wise HAP/BIM/FAM orchestration. gbits owns HAP v1,
# gmat owns validated BIM/FAM metadata, and gsim owns alignment/publication.

.gsim_hap_targets <- function(prefix) {
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !nzchar(prefix)) {
    .gsim_stop("prefix must be one nonempty path string.")
  }
  parent <- normalizePath(dirname(prefix), winslash = "/", mustWork = TRUE)
  base <- basename(prefix)
  if (!nzchar(base) || grepl("\\.(hap|bim|fam)$", base, ignore.case = TRUE)) {
    .gsim_stop("prefix must be an extension-free phased-dataset prefix.")
  }
  canonical <- file.path(parent, base)
  stats::setNames(paste0(canonical, c(".hap", ".bim", ".fam")),
                  c("hap", "bim", "fam"))
}

.gsim_hap_sink_create <- function(backend, path, sample_ids,
                                  overwrite = FALSE, provenance = list()) {
  if (!inherits(backend, "gsim_gbits_backend")) {
    .gsim_stop("backend must be created by .gsim_gbits_backend().")
  }
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    .gsim_stop("path must be one nonempty string.")
  }
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  path <- file.path(parent, basename(path))
  sample_ids <- enc2utf8(as.character(sample_ids))
  if (!length(sample_ids) || anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    .gsim_stop("sample_ids must be nonempty, unique, nonmissing identifiers.")
  }
  overwrite <- .gsim_bed_sink_flag(overwrite, "overwrite")
  if (!is.list(provenance)) .gsim_stop("provenance must be a list.")
  state <- new.env(parent = emptyenv())
  state$pointer <- .Call(C_gsim_gbits_hap_sink_create, backend, enc2utf8(path),
                         length(sample_ids), overwrite)
  state$backend <- backend
  state$path <- path
  state$sample_ids <- sample_ids
  state$chromosome <- character()
  state$marker_count <- integer()
  state$variant_ids <- character()
  state$provenance <- provenance
  state$finalized <- FALSE
  state$cancelled <- FALSE
  class(state) <- "gsim_hap_sink"
  state
}

.gsim_hap_sink_append <- function(sink, chromosome, h1, h2,
                                  variant_ids = attr(h1, "variant_ids", exact = TRUE)) {
  if (!inherits(sink, "gsim_hap_sink") || !is.environment(sink) ||
      is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental HAP sink.")
  }
  if (sink$finalized || sink$cancelled) {
    .gsim_stop("Cannot append after HAP sink finalization or cancellation.")
  }
  chromosome <- enc2utf8(as.character(chromosome))
  if (length(chromosome) != 1L || is.na(chromosome) || !nzchar(chromosome) ||
      chromosome %in% sink$chromosome) {
    .gsim_stop("chromosome must be one new nonempty exact label.")
  }
  if (!inherits(h1, "gsim_gbits_haplotypes") ||
      !inherits(h2, "gsim_gbits_haplotypes")) {
    .gsim_stop("h1 and h2 must be packed gbits haplotype handles.")
  }
  info1 <- .gsim_gbits_info(h1)
  info2 <- .gsim_gbits_info(h2)
  if (!identical(unname(info1[1:2]), unname(info2[1:2]))) {
    .gsim_stop("Packed H1 and H2 dimensions must be identical.")
  }
  if (!identical(attr(h1, "sample_ids", exact = TRUE), sink$sample_ids) ||
      !identical(attr(h2, "sample_ids", exact = TRUE), sink$sample_ids)) {
    .gsim_stop("Packed H1/H2 sample order must exactly match HAP sample_ids.")
  }
  variant_ids <- enc2utf8(as.character(variant_ids))
  marker_count <- as.integer(info1[[2L]])
  if (length(variant_ids) != marker_count || anyNA(variant_ids) ||
      any(!nzchar(variant_ids)) || anyDuplicated(variant_ids) ||
      any(variant_ids %in% sink$variant_ids) ||
      !identical(attr(h1, "variant_ids", exact = TRUE), variant_ids) ||
      !identical(attr(h2, "variant_ids", exact = TRUE), variant_ids)) {
    .gsim_stop("variant_ids must exactly match globally unique H1/H2 order.")
  }
  invisible(.Call(C_gsim_gbits_hap_sink_append, sink$pointer, h1, h2))
  sink$chromosome <- c(sink$chromosome, chromosome)
  sink$marker_count <- c(sink$marker_count, marker_count)
  sink$variant_ids <- c(sink$variant_ids, variant_ids)
  invisible(sink)
}

.gsim_hap_sink_info <- function(sink) {
  if (!inherits(sink, "gsim_hap_sink") || is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental HAP sink.")
  }
  .Call(C_gsim_gbits_hap_sink_info, sink$pointer)
}

.gsim_hap_sink_finalize <- function(sink) {
  if (!inherits(sink, "gsim_hap_sink") || is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental HAP sink.")
  }
  if (sink$finalized) .gsim_stop("HAP sink has already been finalized.")
  if (sink$cancelled) .gsim_stop("A cancelled HAP sink cannot be finalized.")
  invisible(.Call(C_gsim_gbits_hap_sink_finalize, sink$pointer))
  sink$finalized <- TRUE
  native <- .gsim_hap_sink_info(sink)
  words <- ceiling(length(sink$sample_ids) / 64)
  payload <- 2 * sum(sink$marker_count) * words * 8
  expected <- 64 + payload + 48 * length(sink$chromosome)
  list(path = sink$path, individual_count = length(sink$sample_ids),
       marker_count = length(sink$variant_ids),
       chromosome_count = length(sink$chromosome),
       bytes_written = unname(native[["bytes_written"]]),
       expected_bytes = expected, header_bytes = 64,
       chromosome_table_bytes = 48 * length(sink$chromosome),
       packed_data_bytes = payload, state = attr(native, "state", exact = TRUE))
}

.gsim_hap_sink_cancel <- function(sink) {
  if (!inherits(sink, "gsim_hap_sink") || is.null(sink$pointer)) {
    .gsim_stop("sink is not a valid experimental HAP sink.")
  }
  if (sink$finalized) .gsim_stop("A finalized HAP sink cannot be cancelled.")
  if (sink$cancelled) return(invisible(sink))
  invisible(.Call(C_gsim_gbits_hap_sink_cancel, sink$pointer))
  sink$cancelled <- TRUE
  invisible(sink)
}

.gsim_hap_dataset_create <- function(backend, metadata_backend, prefix,
                                     sample_metadata, overwrite = FALSE,
                                     provenance = list()) {
  samples <- .gsim_plink_validate_samples(metadata_backend, sample_metadata)
  targets <- .gsim_hap_targets(prefix)
  overwrite <- .gsim_bed_sink_flag(overwrite, "overwrite")
  exists <- file.exists(targets)
  if (!overwrite && any(exists)) {
    .gsim_stop("Phased dataset destination exists and overwrite is disabled: ",
               paste(targets[exists], collapse = ", "))
  }
  if (overwrite && any(exists) && !all(exists)) {
    .gsim_stop("Overwrite requires an existing complete HAP/BIM/FAM triplet.")
  }
  stage_prefix <- tempfile(pattern = ".gsim-hap-stage-",
                           tmpdir = dirname(targets[[1L]]))
  staged <- stats::setNames(paste0(stage_prefix, c(".hap", ".bim", ".fam")),
                            names(targets))
  state <- new.env(parent = emptyenv())
  state$hap <- .gsim_hap_sink_create(
    backend, staged[["hap"]], samples$metadata$individual_id,
    overwrite = FALSE, provenance = provenance
  )
  state$backend <- backend
  state$metadata_backend <- metadata_backend
  state$sample_metadata <- samples$metadata
  state$sample_pointer <- samples$pointer
  state$targets <- targets
  state$staged <- staged
  state$overwrite <- overwrite
  state$variant_metadata <- list()
  state$provenance <- provenance
  state$finalized <- FALSE
  state$cancelled <- FALSE
  state$failed <- FALSE
  class(state) <- "gsim_hap_dataset"
  state
}

.gsim_hap_dataset_append <- function(dataset, chromosome, h1, h2,
                                     variant_metadata) {
  if (!inherits(dataset, "gsim_hap_dataset") || !is.environment(dataset) ||
      is.null(dataset$hap)) {
    .gsim_stop("dataset is not a valid experimental phased dataset sink.")
  }
  if (dataset$finalized || dataset$cancelled || dataset$failed) {
    .gsim_stop("Cannot append after dataset finalization, cancellation, or failure.")
  }
  metadata <- .gsim_plink_normalize_variants(chromosome, variant_metadata)
  if (!identical(metadata$variant_id, attr(h1, "variant_ids", exact = TRUE)) ||
      !identical(metadata$variant_id, attr(h2, "variant_ids", exact = TRUE))) {
    .gsim_stop("Variant metadata order must exactly match packed H1/H2 IDs.")
  }
  invisible(.gsim_gmat_variant_pointer(dataset$metadata_backend, metadata))
  .gsim_hap_sink_append(dataset$hap, chromosome, h1, h2, metadata$variant_id)
  dataset$variant_metadata[[length(dataset$variant_metadata) + 1L]] <- metadata
  invisible(dataset)
}

.gsim_hap_dataset_finalize <- function(dataset, .test_fail_stage = "none",
                                       .test_fail_publish_after = 0L) {
  if (!inherits(dataset, "gsim_hap_dataset") || !is.environment(dataset)) {
    .gsim_stop("dataset is not a valid experimental phased dataset sink.")
  }
  if (dataset$finalized) .gsim_stop("Phased dataset has already been finalized.")
  if (dataset$cancelled) .gsim_stop("A cancelled phased dataset cannot be finalized.")
  if (dataset$failed) .gsim_stop("A failed phased dataset cannot be finalized.")
  if (!length(dataset$variant_metadata)) {
    .gsim_stop("At least one chromosome must be appended before finalization.")
  }
  if (!.test_fail_stage %in% c("none", "after_hap", "after_bim")) {
    .gsim_stop(".test_fail_stage must be none, after_hap, or after_bim.")
  }
  success <- FALSE
  on.exit({
    if (!success) {
      dataset$failed <- TRUE
      if (!dataset$hap$finalized && !dataset$hap$cancelled) {
        try(.gsim_hap_sink_cancel(dataset$hap), silent = TRUE)
      }
      unlink(c(unname(dataset$staged), paste0(unname(dataset$staged), ".backup")),
             force = TRUE)
    }
  }, add = TRUE)
  variants <- do.call(rbind, dataset$variant_metadata)
  rownames(variants) <- NULL
  variant_pointer <- .gsim_gmat_variant_pointer(dataset$metadata_backend, variants)
  hap_manifest <- .gsim_hap_sink_finalize(dataset$hap)
  if (.test_fail_stage == "after_hap") .gsim_stop("injected failure after HAP completion")
  bim_info <- .Call(C_gsim_gmat_write_bim, variant_pointer,
                    enc2utf8(dataset$staged[["bim"]]))
  if (.test_fail_stage == "after_bim") .gsim_stop("injected failure after BIM completion")
  fam_info <- .Call(C_gsim_gmat_write_fam, dataset$sample_pointer,
                    enc2utf8(dataset$staged[["fam"]]))
  observed <- unname(file.info(dataset$staged)$size)
  if (anyNA(observed) || observed[[1L]] != hap_manifest$expected_bytes ||
      unname(bim_info[["record_count"]]) != nrow(variants) ||
      unname(fam_info[["record_count"]]) != nrow(dataset$sample_metadata)) {
    .gsim_stop("Staged HAP/BIM/FAM triplet failed count or byte validation.")
  }
  .gsim_dataset_publish(dataset$staged, dataset$targets, dataset$overwrite,
                        .test_fail_after = .test_fail_publish_after)
  dataset$finalized <- TRUE
  success <- TRUE
  marker_count <- vapply(dataset$variant_metadata, nrow, integer(1))
  chromosome <- vapply(dataset$variant_metadata, function(x) x$chromosome[[1L]],
                       character(1))
  start <- cumsum(c(1L, head(marker_count, -1L)))
  manifest <- list(
    paths = dataset$targets, individual_count = nrow(dataset$sample_metadata),
    marker_count = nrow(variants), format = "HAP v1", signature = "HAP\\x01",
    chromosomes = data.frame(
      chromosome = chromosome, marker_count = marker_count,
      start_marker = start, end_marker = start + marker_count - 1L,
      stringsAsFactors = FALSE),
    expected_hap_bytes = hap_manifest$expected_bytes,
    observed_hap_bytes = unname(file.info(dataset$targets[["hap"]])$size),
    header_bytes = hap_manifest$header_bytes,
    chromosome_table_bytes = hap_manifest$chromosome_table_bytes,
    packed_data_bytes = hap_manifest$packed_data_bytes,
    bim_record_count = unname(bim_info[["record_count"]]),
    fam_record_count = unname(fam_info[["record_count"]]),
    sample_ids = dataset$sample_metadata$individual_id,
    variant_ids = variants$variant_id,
    allele_orientation = "bit 1 = ALT = BIM A1; bit 0 = REF = BIM A2",
    backend = list(
      gbits = list(version = attr(dataset$backend, "gbits_version", exact = TRUE),
                   abi = attr(dataset$backend, "gbits_abi", exact = TRUE)),
      gmat = list(version = attr(dataset$metadata_backend, "gmat_version", exact = TRUE),
                  abi = attr(dataset$metadata_backend, "gmat_abi", exact = TRUE))),
    provenance = dataset$provenance, publication_status = "published",
    transaction = "same-directory staging with backup-and-rollback publication"
  )
  class(manifest) <- c("gsim_hap_dataset_manifest", "list")
  manifest
}

.gsim_hap_dataset_cancel <- function(dataset) {
  if (!inherits(dataset, "gsim_hap_dataset") || !is.environment(dataset)) {
    .gsim_stop("dataset is not a valid experimental phased dataset sink.")
  }
  if (dataset$finalized) .gsim_stop("A finalized phased dataset cannot be cancelled.")
  if (dataset$cancelled) return(invisible(dataset))
  if (!dataset$hap$finalized && !dataset$hap$cancelled) {
    .gsim_hap_sink_cancel(dataset$hap)
  }
  unlink(c(unname(dataset$staged), paste0(unname(dataset$staged), ".backup")),
         force = TRUE)
  dataset$cancelled <- TRUE
  invisible(dataset)
}

.gsim_hap_dataset_open <- function(backend, metadata_backend, prefix) {
  if (!inherits(backend, "gsim_gbits_backend") ||
      !inherits(metadata_backend, "gsim_gmat_backend")) {
    .gsim_stop("backend and metadata_backend must be gbits/gmat backends.")
  }
  targets <- .gsim_hap_targets(prefix)
  if (!all(file.exists(targets))) {
    .gsim_stop("A complete HAP/BIM/FAM triplet is required.")
  }
  bim <- .Call(C_gsim_gmat_read_bim, metadata_backend, enc2utf8(targets[["bim"]]))
  fam <- .Call(C_gsim_gmat_read_fam, metadata_backend, enc2utf8(targets[["fam"]]))
  variants <- data.frame(
    chromosome = bim$chromosome, variant_id = bim$variant_id,
    genetic_position_cm = bim$genetic_position_cm,
    base_pair_position = bim$base_pair_position, alt = bim$alt, ref = bim$ref,
    bit1_allele = bim$alt, bit0_allele = bim$ref, stringsAsFactors = FALSE)
  samples <- data.frame(
    family_id = fam$family_id, individual_id = fam$individual_id,
    paternal_id = fam$paternal_id, maternal_id = fam$maternal_id,
    sex = fam$sex, phenotype = rep.int(NA_real_, length(fam$individual_id)),
    stringsAsFactors = FALSE)
  pointer <- .Call(C_gsim_gbits_hap_reader_open, backend, enc2utf8(targets[["hap"]]))
  info <- .Call(C_gsim_gbits_hap_reader_info, pointer)
  runs <- rle(variants$chromosome)
  expected_start <- cumsum(c(0, head(runs$lengths, -1L)))
  ranges <- info$ranges
  if (info$individual_count != nrow(samples) || info$marker_count != nrow(variants) ||
      info$chromosome_count != length(runs$values) ||
      !identical(as.double(ranges[, "global_start_marker"]), as.double(expected_start)) ||
      !identical(as.double(ranges[, "marker_count"]), as.double(runs$lengths))) {
    try(.Call(C_gsim_gbits_hap_reader_close, pointer), silent = TRUE)
    .gsim_stop("HAP dimensions/ranges do not align with BIM/FAM metadata.")
  }
  state <- new.env(parent = emptyenv())
  state$pointer <- pointer
  state$backend <- backend
  state$metadata_backend <- metadata_backend
  state$paths <- targets
  state$info <- info
  state$variants <- variants
  state$samples <- samples
  state$chromosome <- runs$values
  state$marker_count <- runs$lengths
  state$closed <- FALSE
  class(state) <- "gsim_hap_dataset_reader"
  state
}

.gsim_hap_dataset_inspect <- function(dataset) {
  if (!inherits(dataset, "gsim_hap_dataset_reader") || dataset$closed) {
    .gsim_stop("dataset is not an open phased dataset.")
  }
  start <- cumsum(c(1L, head(dataset$marker_count, -1L)))
  list(paths = dataset$paths,
       individual_count = nrow(dataset$samples),
       marker_count = nrow(dataset$variants),
       chromosome_count = length(dataset$chromosome),
       format = "HAP v1", signature = "HAP\\x01",
       chromosomes = data.frame(
         chromosome = dataset$chromosome, marker_count = dataset$marker_count,
         start_marker = start, end_marker = start + dataset$marker_count - 1L,
         stringsAsFactors = FALSE),
       sample_ids = dataset$samples$individual_id,
       variant_ids = dataset$variants$variant_id,
       allele_orientation = "bit 1 = ALT = BIM A1; bit 0 = REF = BIM A2",
       validation = "counts and contiguous BIM blocks; v1 has no metadata checksum")
}

.gsim_hap_dataset_load_chromosome <- function(dataset, chromosome) {
  if (!inherits(dataset, "gsim_hap_dataset_reader") || dataset$closed) {
    .gsim_stop("dataset is not an open phased dataset.")
  }
  chromosome <- enc2utf8(as.character(chromosome))
  index <- match(chromosome, dataset$chromosome)
  if (length(chromosome) != 1L || is.na(chromosome) || is.na(index)) {
    .gsim_stop("chromosome must exactly identify one HAP chromosome.")
  }
  loaded <- .Call(C_gsim_gbits_hap_reader_load, dataset$pointer, as.integer(index))
  start <- sum(head(dataset$marker_count, index - 1L)) + 1L
  marker_ids <- dataset$variants$variant_id[
    seq.int(start, length.out = dataset$marker_count[[index]])]
  loaded$h1 <- .gsim_gbits_tag(loaded$h1, dataset$samples$individual_id, marker_ids)
  loaded$h2 <- .gsim_gbits_tag(loaded$h2, dataset$samples$individual_id, marker_ids)
  attr(loaded, "chromosome") <- chromosome
  loaded
}

.gsim_hap_dataset_close <- function(dataset) {
  if (!inherits(dataset, "gsim_hap_dataset_reader")) {
    .gsim_stop("dataset is not a phased dataset reader.")
  }
  if (dataset$closed) return(invisible(dataset))
  invisible(.Call(C_gsim_gbits_hap_reader_close, dataset$pointer))
  dataset$closed <- TRUE
  invisible(dataset)
}
