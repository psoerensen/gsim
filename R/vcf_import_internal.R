# Strict phased-VCF import: gsim owns parsing, identity, and packed bits,
# and gsim owns map/sample alignment and HAP/BIM/FAM publication.

.gsim_vcf_reader_open <- function(metadata_backend, vcf) {
  if (!inherits(metadata_backend, "gsim_metadata_backend")) {
    .gsim_stop("metadata_backend must be created by .gsim_metadata_backend().")
  }
  if (!is.character(vcf) || length(vcf) != 1L || is.na(vcf) || !nzchar(vcf)) {
    .gsim_stop("vcf must be one nonempty path string.")
  }
  path <- normalizePath(vcf, winslash = "/", mustWork = TRUE)
  value <- .Call(C_gsim_metadata_vcf_open, metadata_backend, enc2utf8(path))
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
  if (nrow(map) != nrow(variants) || anyNA(chromosome) ||
      any(!nzchar(chromosome))) {
    .gsim_stop("map must contain one nonmissing row per imported variant.")
  }
  if (has_id) {
    key <- enc2utf8(as.character(map$variant_id))
    imported_key <- variants$variant_id
  } else {
    bp <- suppressWarnings(as.double(map$base_pair_position))
    if (anyNA(bp) || any(!is.finite(bp)) || any(bp <= 0) ||
        any(bp != floor(bp)) || any(bp > 2^53 - 1)) {
      .gsim_stop("map base_pair_position must contain exact positive integers.")
    }
    key <- paste(chromosome, format(bp, scientific = FALSE, trim = TRUE), sep = "\r")
    imported_key <- paste(variants$chromosome,
                          format(variants$base_pair_position,
                                 scientific = FALSE, trim = TRUE), sep = "\r")
  }
  if (anyNA(key) || any(!nzchar(key)) || anyDuplicated(key) ||
      !setequal(key, imported_key)) {
    .gsim_stop("map alignment keys must match every imported variant exactly once.")
  }
  aligned <- map[match(imported_key, key), , drop = FALSE]
  if (!identical(enc2utf8(as.character(aligned$chromosome)),
                 variants$chromosome)) {
    .gsim_stop("map chromosome labels do not match imported variants exactly.")
  }
  cm <- suppressWarnings(as.double(aligned$genetic_position_cm))
  if (anyNA(cm) || any(!is.finite(cm)) || any(cm < 0)) {
    .gsim_stop("map genetic_position_cm must be finite and nonnegative.")
  }
  runs <- rle(variants$chromosome)
  start <- cumsum(c(1L, head(runs$lengths, -1L)))
  for (i in seq_along(start)) {
    rows <- seq.int(start[[i]], length.out = runs$lengths[[i]])
    if (any(diff(cm[rows]) < 0)) {
      .gsim_stop("map genetic_position_cm must be nondecreasing within chromosome.")
    }
  }
  cm
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
  overwrite = FALSE
) {
  reader <- .gsim_vcf_reader_open(metadata_backend, vcf)
  on.exit(try(.gsim_vcf_reader_close(reader), silent = TRUE), add = TRUE)
  variants <- reader$variants
  variants$genetic_position_cm <- .gsim_vcf_align_map(map, variants)
  samples <- .gsim_vcf_sample_metadata(reader$samples, sample_metadata)
  provenance <- list(
    operation = "strict phased biallelic VCF import",
    source_vcf = reader$path,
    vcf_subset = "plain text; diploid phased GT; biallelic uppercase A/C/G/T SNPs",
    generated_variant_ids = variants$variant_id[variants$generated_id],
    map_alignment = if ("variant_id" %in% names(map)) "exact variant_id" else
      "exact chromosome and base_pair_position",
    allele_orientation = "VCF REF = bit 0 = BIM A2; VCF ALT = bit 1 = BIM A1",
    phase_orientation = "GT left allele = H1; GT right allele = H2"
  )
  dataset <- .gsim_hap_dataset_create(
    backend, metadata_backend, output, samples, overwrite, provenance
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
  manifest$import <- list(
    vcf_bytes = unname(file.info(reader$path)$size),
    peak_chromosome_packed_bytes = peak_packed,
    native_record_allele_buffer_bytes = 2 * length(reader$samples),
    r_record_allele_buffer_bytes = 2 * length(reader$samples),
    maximum_record_allele_buffer_bytes = 4 * length(reader$samples),
    dense_haplotype_bytes_avoided = 2 * nrow(variants) * length(reader$samples),
    dense_genotype_bytes_avoided = nrow(variants) * length(reader$samples),
    generated_variant_ids = variants$variant_id[variants$generated_id]
  )
  manifest
}
