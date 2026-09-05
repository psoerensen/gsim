# Strict phased-VCF import: gsim owns parsing, identity, and packed bits,
# and gsim owns map/sample alignment and HAP/BIM/FAM publication.

.gsim_vcf_reader_open <- function(
  metadata_backend, vcf, samples = NULL, chromosome = NULL, region = NULL,
  unsupported = "error"
) {
  if (!inherits(metadata_backend, "gsim_metadata_backend")) {
    .gsim_stop("metadata_backend must be created by .gsim_metadata_backend().")
  }
  if (!is.character(vcf) || length(vcf) != 1L || is.na(vcf) || !nzchar(vcf)) {
    .gsim_stop("vcf must be one nonempty path string.")
  }
  if (!is.null(samples)) {
    samples <- enc2utf8(as.character(samples))
    if (!length(samples) || anyNA(samples) || any(!nzchar(samples)) ||
        anyDuplicated(samples)) {
      .gsim_stop("samples must contain unique, nonmissing, nonempty VCF sample IDs.")
    }
  }
  if (!is.null(chromosome)) {
    chromosome <- enc2utf8(as.character(chromosome))
    if (length(chromosome) != 1L || is.na(chromosome) || !nzchar(chromosome)) {
      .gsim_stop("chromosome must be NULL or one nonempty exact label.")
    }
  }
  if (!is.null(region)) {
    region <- suppressWarnings(as.double(region))
    if (length(region) != 2L || anyNA(region) || any(!is.finite(region)) ||
        any(region <= 0) || any(region != floor(region)) ||
        any(region > 2^53 - 1) || region[[1L]] > region[[2L]] ||
        is.null(chromosome)) {
      .gsim_stop("region requires chromosome and two exact positive integers with start <= end.")
    }
  }
  unsupported <- match.arg(unsupported, c("skip", "error"))
  path <- normalizePath(vcf, winslash = "/", mustWork = TRUE)
  value <- .Call(
    C_gsim_metadata_vcf_open, metadata_backend, enc2utf8(path), samples,
    chromosome, region, unsupported
  )
  value$variants <- as.data.frame(value$variants, stringsAsFactors = FALSE)
  value$chromosomes <- as.data.frame(value$chromosomes, stringsAsFactors = FALSE)
  value$path <- path
  value$closed <- FALSE
  class(value) <- c("gsim_vcf_reader", "list")
  value
}

.gsim_vcf_reader_close <- function(reader) {
  if (!inherits(reader, "gsim_vcf_reader")) {
    .gsim_stop("reader is not a gsim VCF reader.")
  }
  if (reader$closed) return(invisible(reader))
  invisible(.Call(C_gsim_metadata_vcf_close, reader$pointer))
  reader$closed <- TRUE
  invisible(reader)
}

.gsim_vcf_reader_start <- function(reader, chromosome_index) {
  if (!inherits(reader, "gsim_vcf_reader") || reader$closed) {
    .gsim_stop("reader is not an open gsim VCF reader.")
  }
  invisible(.Call(C_gsim_metadata_vcf_start, reader$pointer,
                  as.integer(chromosome_index)))
}

.gsim_vcf_reader_next <- function(reader) {
  if (!inherits(reader, "gsim_vcf_reader") || reader$closed) {
    .gsim_stop("reader is not an open gsim VCF reader.")
  }
  .Call(C_gsim_metadata_vcf_next, reader$pointer)
}

.gsim_vcf_align_map <- function(map, variants) {
  if (!is.data.frame(map)) .gsim_stop("map must be a data frame.")
  required <- c("chromosome", "genetic_position_cm")
  if (any(!required %in% names(map))) {
    .gsim_stop("map must contain chromosome and genetic_position_cm.")
  }
  has_id <- "variant_id" %in% names(map)
  has_bp <- "base_pair_position" %in% names(map)
  if (has_id == has_bp) {
    .gsim_stop("map must contain exactly one alignment key: variant_id or base_pair_position.")
  }
  chromosome <- enc2utf8(as.character(map$chromosome))
  if (!nrow(map) || anyNA(chromosome) || any(!nzchar(chromosome))) {
    .gsim_stop("map chromosome labels must be nonmissing and nonempty.")
  }
  cm <- suppressWarnings(as.double(map$genetic_position_cm))
  if (anyNA(cm) || any(!is.finite(cm)) || any(cm < 0)) {
    .gsim_stop("map genetic_position_cm must be finite and nonnegative.")
  }
  if (has_id) {
    key <- enc2utf8(as.character(map$variant_id))
    imported_key <- variants$variant_id
    if (nrow(map) != nrow(variants) || anyNA(key) || any(!nzchar(key)) ||
        anyDuplicated(key) || !setequal(key, imported_key)) {
      .gsim_stop("map variant_id alignment must match every imported variant exactly once.")
    }
    aligned <- map[match(imported_key, key), , drop = FALSE]
    if (!identical(enc2utf8(as.character(aligned$chromosome)),
                   variants$chromosome)) {
      .gsim_stop("map chromosome labels do not match imported variants exactly.")
    }
    aligned_cm <- suppressWarnings(as.double(aligned$genetic_position_cm))
    runs <- rle(variants$chromosome)
    start <- cumsum(c(1L, head(runs$lengths, -1L)))
    for (i in seq_along(start)) {
      rows <- seq.int(start[[i]], length.out = runs$lengths[[i]])
      if (any(diff(aligned_cm[rows]) < 0)) {
        .gsim_stop("map genetic_position_cm must be nondecreasing within chromosome.")
      }
    }
    attr(aligned_cm, "interpolation_policy") <- "exact variant_id alignment; no interpolation"
    return(aligned_cm)
  }

  bp <- suppressWarnings(as.double(map$base_pair_position))
  if (anyNA(bp) || any(!is.finite(bp)) || any(bp <= 0) ||
      any(bp != floor(bp)) || any(bp > 2^53 - 1)) {
    .gsim_stop("map base_pair_position must contain exact positive integers.")
  }
  map_runs <- rle(chromosome)
  if (anyDuplicated(map_runs$values)) {
    .gsim_stop("map chromosome labels must occupy one contiguous block.")
  }
  imported_chromosomes <- unique(variants$chromosome)
  if (!setequal(map_runs$values, imported_chromosomes)) {
    .gsim_stop("sparse map chromosome labels must match imported chromosomes exactly.")
  }
  map_start <- cumsum(c(1L, head(map_runs$lengths, -1L)))
  for (i in seq_along(map_start)) {
    rows <- seq.int(map_start[[i]], length.out = map_runs$lengths[[i]])
    if (any(diff(bp[rows]) <= 0)) {
      .gsim_stop("map base_pair_position must be strictly increasing within chromosome.")
    }
    if (any(diff(cm[rows]) < 0)) {
      .gsim_stop("map genetic_position_cm must be nondecreasing within chromosome.")
    }
  }
  interpolated <- numeric(nrow(variants))
  for (label in unique(variants$chromosome)) {
    target_rows <- which(variants$chromosome == label)
    knot_rows <- which(chromosome == label)
    knot_bp <- bp[knot_rows]
    knot_cm <- cm[knot_rows]
    target_bp <- variants$base_pair_position[target_rows]
    if (length(knot_rows) == 1L && all(target_bp == knot_bp)) {
      interpolated[target_rows] <- knot_cm
      next
    }
    if (length(knot_rows) < 2L) {
      .gsim_stop(paste0("sparse map chromosome '", label,
                        "' requires two knots unless every retained variant is at its sole knot."))
    }
    if (any(target_bp < knot_bp[[1L]]) ||
        any(target_bp > utils::tail(knot_bp, 1L))) {
      .gsim_stop(paste0("retained variants on chromosome '", label,
                        "' fall outside the supplied map range; extrapolation is not allowed."))
    }
    left <- findInterval(target_bp, knot_bp)
    exact <- knot_bp[left] == target_bp
    values <- knot_cm[left]
    between <- which(!exact)
    if (length(between)) {
      lo <- left[between]
      hi <- lo + 1L
      fraction <- (target_bp[between] - knot_bp[lo]) /
        (knot_bp[hi] - knot_bp[lo])
      values[between] <- knot_cm[lo] +
        (knot_cm[hi] - knot_cm[lo]) * fraction
    }
    interpolated[target_rows] <- values
  }
  attr(interpolated, "interpolation_policy") <- paste(
    "piecewise linear cumulative-cM interpolation:",
    "exact knots copied; otherwise cm_left + (cm_right - cm_left) *",
    "(bp - bp_left) / (bp_right - bp_left); no extrapolation"
  )
  interpolated
}

.gsim_vcf_sample_metadata <- function(samples, sample_metadata = NULL) {
  samples <- enc2utf8(as.character(samples))
  if (is.null(sample_metadata)) {
    return(.gsim_plink_sample_metadata(
      samples, family_id = rep.int("reference", length(samples))
    ))
  }
  if (!is.data.frame(sample_metadata) ||
      !"individual_id" %in% names(sample_metadata)) {
    .gsim_stop("sample_metadata must be a data frame containing individual_id.")
  }
  ids <- enc2utf8(as.character(sample_metadata$individual_id))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      !setequal(ids, samples) || nrow(sample_metadata) != length(samples)) {
    .gsim_stop("sample_metadata individual_id must match every VCF sample exactly once.")
  }
  value <- sample_metadata[match(samples, ids), , drop = FALSE]
  field <- function(name, default) {
    if (name %in% names(value)) value[[name]] else rep.int(default, nrow(value))
  }
  out <- .gsim_plink_sample_metadata(
    samples,
    family_id = field("family_id", "reference"),
    paternal_id = field("paternal_id", "0"),
    maternal_id = field("maternal_id", "0"),
    sex = field("sex", 0L),
    phenotype = field("phenotype", NA_real_)
  )
  if (any(as.character(out$paternal_id) != "0") ||
      any(as.character(out$maternal_id) != "0")) {
    .gsim_stop("Imported reference samples must have paternal_id and maternal_id 0.")
  }
  out
}

.gsim_import_vcf_internal <- function(
  backend, metadata_backend, vcf, map, output, sample_metadata = NULL,
  samples = NULL, chromosome = NULL, region = NULL,
  unsupported = "error", overwrite = FALSE
) {
  reader <- .gsim_vcf_reader_open(
    metadata_backend, vcf, samples, chromosome, region, unsupported
  )
  on.exit(try(.gsim_vcf_reader_close(reader), silent = TRUE), add = TRUE)
  variants <- reader$variants
  aligned_map <- .gsim_vcf_align_map(map, variants)
  interpolation_policy <- attr(aligned_map, "interpolation_policy")
  variants$genetic_position_cm <- unname(aligned_map)
  sample_table <- .gsim_vcf_sample_metadata(reader$samples, sample_metadata)
  provenance <- list(
    operation = "streaming phased biallelic VCF import",
    source_vcf = reader$path,
    input_type = reader$report$input_type,
    vcf_subset = paste(
      "plain/gzip/BGZF VCF; selected diploid phased GT;",
      "biallelic uppercase A/C/G/T SNPs"
    ),
    unsupported = unsupported,
    selected_chromosome = chromosome,
    selected_region = region,
    generated_variant_id_count = sum(variants$generated_id),
    map_alignment = interpolation_policy,
    allele_orientation = "VCF REF = bit 0 = BIM A2; VCF ALT = bit 1 = BIM A1",
    phase_orientation = "GT left allele = H1; GT right allele = H2"
  )
  dataset <- .gsim_hap_dataset_create(
    backend, metadata_backend, output, sample_table, overwrite, provenance
  )
  completed <- FALSE
  on.exit({
    if (!completed) try(.gsim_hap_dataset_cancel(dataset), silent = TRUE)
  }, add = TRUE)
  peak_packed <- 0
  for (block_index in seq_len(nrow(reader$chromosomes))) {
    block <- reader$chromosomes[block_index, , drop = FALSE]
    first <- as.integer(block$first_variant)
    count <- as.integer(block$variant_count)
    rows <- seq.int(first, length.out = count)
    marker_ids <- variants$variant_id[rows]
    h1 <- .gsim_packed_zero(backend, length(reader$samples), count,
                           reader$samples, marker_ids)
    h2 <- .gsim_packed_zero(backend, length(reader$samples), count,
                           reader$samples, marker_ids)
    on.exit({
      try(.gsim_packed_close(h1), silent = TRUE)
      try(.gsim_packed_close(h2), silent = TRUE)
    }, add = TRUE)
    .gsim_vcf_reader_start(reader, block_index)
    for (local_marker in seq_len(count)) {
      record <- .gsim_vcf_reader_next(reader)
      if (is.null(record) || record$variant_index != rows[[local_marker]]) {
        .gsim_stop("VCF chromosome stream did not preserve marker order.")
      }
      .gsim_packed_set_marker(h1, h2, local_marker, record$h1, record$h2)
    }
    if (!is.null(.gsim_vcf_reader_next(reader))) {
      .gsim_stop("VCF chromosome stream exceeded its declared block.")
    }
    variant_metadata <- data.frame(
      chromosome = variants$chromosome[rows], variant_id = marker_ids,
      genetic_position_cm = variants$genetic_position_cm[rows],
      base_pair_position = variants$base_pair_position[rows],
      alt = variants$alt[rows], ref = variants$ref[rows],
      bit1_allele = variants$alt[rows], bit0_allele = variants$ref[rows],
      stringsAsFactors = FALSE
    )
    .gsim_hap_dataset_append(dataset, block$chromosome[[1L]], h1, h2,
                             variant_metadata)
    peak_packed <- max(peak_packed,
                       .gsim_packed_info(h1)[["storage_bytes"]] +
                         .gsim_packed_info(h2)[["storage_bytes"]])
    .gsim_packed_close(h1)
    .gsim_packed_close(h2)
    h1 <- h2 <- NULL
  }
  manifest <- .gsim_hap_dataset_finalize(dataset)
  completed <- TRUE
  output_sizes <- unname(file.info(manifest$paths)$size)
  names(output_sizes) <- names(manifest$paths)
  manifest$import <- c(reader$report, list(
    input_path = reader$path,
    input_bytes = unname(file.info(reader$path)$size),
    output_chromosomes = stats::setNames(
      as.integer(reader$chromosomes$variant_count), reader$chromosomes$chromosome
    ),
    output_file_bytes = output_sizes,
    peak_chromosome_packed_bytes = peak_packed,
    native_record_allele_buffer_bytes = 2 * length(reader$samples),
    r_record_allele_buffer_bytes = 2 * length(reader$samples),
    maximum_record_allele_buffer_bytes = 4 * length(reader$samples),
    dense_haplotype_bytes_avoided = 2 * nrow(variants) * length(reader$samples),
    dense_genotype_bytes_avoided = nrow(variants) * length(reader$samples),
    dense_haplotype_matrix_allocated = FALSE,
    dense_genotype_matrix_allocated = FALSE,
    full_uncompressed_temporary_vcf_created = FALSE,
    generated_variant_id_count = sum(variants$generated_id),
    genetic_map_interpolation_policy = interpolation_policy
  ))
  manifest
}
