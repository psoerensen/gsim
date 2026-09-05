# Experimental BED/BIM/FAM dataset orchestration. gbits owns BED coding, gmat
# owns validated metadata serialization, and gsim owns aligned publication.

.gsim_gmat_symbols <- c(
  "gmat_abi_version", "gmat_library_version", "gmat_last_error",
  "gmat_variant_metadata_create", "gmat_variant_metadata_close",
  "gmat_variant_metadata_write_bim", "gmat_variant_metadata_read_bim",
  "gmat_variant_metadata_count", "gmat_variant_metadata_get",
  "gmat_sample_metadata_create", "gmat_sample_metadata_close",
  "gmat_sample_metadata_write_fam", "gmat_sample_metadata_read_fam",
  "gmat_sample_metadata_count", "gmat_sample_metadata_get"
)

.gsim_gmat_backend <- function(
  library = Sys.getenv("GSIM_GMAT_LIBRARY", unset = "")
) {
  if (!is.character(library) || length(library) != 1L || is.na(library) ||
      !nzchar(library)) {
    .gsim_stop(
      "A gmat shared-library path is required via library or GSIM_GMAT_LIBRARY."
    )
  }
  path <- normalizePath(library, winslash = "/", mustWork = TRUE)
  dll <- tryCatch(
    dyn.load(path, local = TRUE, now = TRUE),
    error = function(error) .gsim_stop(
      "Unable to load gmat: ", conditionMessage(error)
    )
  )
  addresses <- tryCatch(
    stats::setNames(lapply(.gsim_gmat_symbols, function(name) {
      getNativeSymbolInfo(name, PACKAGE = dll)$address
    }), .gsim_gmat_symbols),
    error = function(error) .gsim_stop(
      "Incompatible gmat library: ", conditionMessage(error)
    )
  )
  pointer <- .Call(C_gsim_gmat_backend, addresses)
  attr(pointer, "library") <- path
  attr(pointer, "dll") <- dll
  attr(pointer, "gmat_abi") <- 0L
  class(pointer) <- "gsim_gmat_backend"
  pointer
}

.gsim_plink_sample_metadata <- function(
  individual_id,
  family_id = rep.int("gsim", length(individual_id)),
  paternal_id = rep.int("0", length(individual_id)),
  maternal_id = rep.int("0", length(individual_id)),
  sex = rep.int(0L, length(individual_id)),
  phenotype = rep.int(NA_real_, length(individual_id))
) {
  data.frame(
    family_id = family_id, individual_id = individual_id,
    paternal_id = paternal_id, maternal_id = maternal_id,
    sex = sex, phenotype = phenotype, stringsAsFactors = FALSE
  )
}

.gsim_plink_pedigree_metadata <- function(
  pedigree,
  individual_id = pedigree$canonical_order,
  family_id = rep.int("gsim", length(individual_id)),
  sex = rep.int(0L, length(individual_id))
) {
  if (!inherits(pedigree, "gsim_pedigree") ||
      !is.data.frame(pedigree$pedigree) || is.null(pedigree$canonical_order)) {
    .gsim_stop("pedigree must be a gsim_pedigree object.")
  }
  canonical <- as.character(pedigree$canonical_order)
  individual_id <- as.character(individual_id)
  if (!identical(individual_id, canonical)) {
    .gsim_stop("individual_id must exactly equal pedigree canonical_order.")
  }
  tab <- pedigree$pedigree
  tab <- tab[match(canonical, as.character(tab$animal)), , drop = FALSE]
  sire <- as.character(tab$sire)
  dam <- as.character(tab$dam)
  sire[is.na(tab$sire)] <- "0"
  dam[is.na(tab$dam)] <- "0"
  .gsim_plink_sample_metadata(
    individual_id, family_id, sire, dam, sex,
    rep.int(NA_real_, length(individual_id))
  )
}

.gsim_plink_validate_samples <- function(metadata_backend, metadata) {
  if (!inherits(metadata_backend, "gsim_gmat_backend")) {
    .gsim_stop("metadata_backend must be created by .gsim_gmat_backend().")
  }
  if (!is.data.frame(metadata)) {
    .gsim_stop("sample_metadata must be a data frame.")
  }
  required <- c(
    "family_id", "individual_id", "paternal_id", "maternal_id", "sex"
  )
  if (any(!required %in% names(metadata))) {
    .gsim_stop("sample_metadata must contain family_id, individual_id, ",
               "paternal_id, maternal_id, and sex.")
  }
  if (!nrow(metadata)) .gsim_stop("sample_metadata must not be empty.")
  phenotype <- if ("phenotype" %in% names(metadata)) {
    metadata$phenotype
  } else {
    rep.int(NA_real_, nrow(metadata))
  }
  if (length(phenotype) != nrow(metadata) ||
      any(!is.na(phenotype) & phenotype != -9)) {
    .gsim_stop("This milestone supports only missing phenotype metadata (-9/NA).")
  }
  values <- lapply(metadata[required[1:4]], as.character)
  if (any(vapply(values, anyNA, logical(1)))) {
    .gsim_stop("FAM identity and parent fields must not be missing; use '0' parents.")
  }
  sex <- suppressWarnings(as.integer(metadata$sex))
  if (length(sex) != nrow(metadata) || anyNA(sex) ||
      any(sex < 0L | sex > 2L) ||
      any(as.character(sex) != as.character(metadata$sex))) {
    .gsim_stop("sex must contain only integer codes 0, 1, or 2.")
  }
  normalized <- data.frame(
    family_id = enc2utf8(values[[1L]]),
    individual_id = enc2utf8(values[[2L]]),
    paternal_id = enc2utf8(values[[3L]]),
    maternal_id = enc2utf8(values[[4L]]),
    sex = sex, phenotype = rep.int(NA_real_, nrow(metadata)),
    stringsAsFactors = FALSE
  )
  pointer <- .Call(
    C_gsim_gmat_sample_create, metadata_backend, normalized$family_id,
    normalized$individual_id, normalized$paternal_id,
    normalized$maternal_id, normalized$sex
  )
  list(metadata = normalized, pointer = pointer)
}

.gsim_plink_normalize_variants <- function(chromosome, metadata) {
  chromosome <- as.character(chromosome)
  if (length(chromosome) != 1L || is.na(chromosome) || !nzchar(chromosome)) {
    .gsim_stop("chromosome must be one nonempty exact label.")
  }
  if (!is.data.frame(metadata)) {
    .gsim_stop("variant_metadata must be a data frame.")
  }
  required <- c(
    "variant_id", "genetic_position_cm", "base_pair_position", "alt", "ref",
    "bit1_allele", "bit0_allele"
  )
  if (any(!required %in% names(metadata))) {
    .gsim_stop("variant_metadata must contain variant_id, genetic_position_cm, ",
               "base_pair_position, alt, ref, bit1_allele, and bit0_allele.")
  }
  if (!nrow(metadata)) .gsim_stop("variant_metadata must not be empty.")
  if ("chromosome" %in% names(metadata)) {
    supplied <- as.character(metadata$chromosome)
    if (anyNA(supplied) || any(supplied != chromosome)) {
      .gsim_stop("variant_metadata chromosome must exactly match chromosome.")
    }
  }
  character_columns <- lapply(metadata[c("variant_id", "alt", "ref",
                                         "bit1_allele", "bit0_allele")],
                              as.character)
  if (any(vapply(character_columns, anyNA, logical(1)))) {
    .gsim_stop("Variant IDs and allele fields must not be missing.")
  }
  names(character_columns) <- c("variant_id", "alt", "ref",
                                "bit1_allele", "bit0_allele")
  if (!identical(character_columns$bit1_allele, character_columns$alt) ||
      !identical(character_columns$bit0_allele, character_columns$ref)) {
    .gsim_stop("Packed allele orientation must explicitly satisfy bit1_allele = ",
               "ALT and bit0_allele = REF; allele reversal is not implicit.")
  }
  genetic <- suppressWarnings(as.double(metadata$genetic_position_cm))
  bp <- suppressWarnings(as.double(metadata$base_pair_position))
  if (length(genetic) != nrow(metadata) || length(bp) != nrow(metadata)) {
    .gsim_stop("Variant metadata columns must have identical lengths.")
  }
  data.frame(
    chromosome = rep.int(enc2utf8(chromosome), nrow(metadata)),
    variant_id = enc2utf8(character_columns$variant_id),
    genetic_position_cm = genetic,
    base_pair_position = bp,
    alt = enc2utf8(character_columns$alt),
    ref = enc2utf8(character_columns$ref),
    bit1_allele = enc2utf8(character_columns$bit1_allele),
    bit0_allele = enc2utf8(character_columns$bit0_allele),
    stringsAsFactors = FALSE
  )
}

.gsim_gmat_variant_pointer <- function(metadata_backend, metadata) {
  .Call(
    C_gsim_gmat_variant_create, metadata_backend, metadata$chromosome,
    metadata$variant_id, metadata$genetic_position_cm,
    metadata$base_pair_position, metadata$alt, metadata$ref
  )
}

.gsim_plink_targets <- function(prefix) {
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !nzchar(prefix)) {
    .gsim_stop("prefix must be one nonempty path string.")
  }
  parent <- normalizePath(dirname(prefix), winslash = "/", mustWork = TRUE)
  base <- basename(prefix)
  if (!nzchar(base) || grepl("\\.(bed|bim|fam)$", base, ignore.case = TRUE)) {
    .gsim_stop("prefix must be an extension-free dataset prefix.")
  }
  canonical <- file.path(parent, base)
  stats::setNames(paste0(canonical, c(".bed", ".bim", ".fam")),
                  c("bed", "bim", "fam"))
}

.gsim_dataset_publish <- function(staged, targets, overwrite,
                                  .test_fail_after = 0L) {
  exists <- file.exists(targets)
  if (!overwrite && any(exists)) {
    .gsim_stop("Dataset destination exists and overwrite is disabled: ",
               paste(targets[exists], collapse = ", "))
  }
  if (overwrite && any(exists) && !all(exists)) {
    .gsim_stop("Overwrite requires an existing complete three-file dataset.")
  }
  if (!all(file.exists(staged))) {
    .gsim_stop("Cannot publish an incomplete staged three-file dataset.")
  }
  fail_after <- .gsim_hapnest_integer_scalar(
    .test_fail_after, ".test_fail_after", 0
  )
  if (fail_after > 3L) .gsim_stop(".test_fail_after must be between 0 and 3.")
  backups <- paste0(staged, ".backup")
  backed_up <- rep.int(FALSE, 3L)
  published <- rep.int(FALSE, 3L)
  rollback <- character()
  result <- tryCatch({
    if (all(exists)) {
      for (i in seq_along(targets)) {
        if (!isTRUE(file.rename(targets[[i]], backups[[i]]))) {
          stop("could not back up ", targets[[i]], call. = FALSE)
        }
        backed_up[[i]] <- TRUE
      }
    }
    for (i in seq_along(targets)) {
      if (!isTRUE(file.rename(staged[[i]], targets[[i]]))) {
        stop("could not publish ", targets[[i]], call. = FALSE)
      }
      published[[i]] <- TRUE
      if (fail_after == i) {
        stop("injected publication failure after final rename ", i,
             call. = FALSE)
      }
    }
    if (any(backed_up)) unlink(backups[backed_up], force = TRUE)
    TRUE
  }, error = function(error) {
    for (i in rev(seq_along(targets))) {
      if (published[[i]] && file.exists(targets[[i]]) &&
          !isTRUE(unlink(targets[[i]], force = TRUE) == 0L)) {
        rollback <<- c(rollback, paste0("remove new ", targets[[i]]))
      }
    }
    for (i in seq_along(targets)) {
      if (backed_up[[i]] && file.exists(backups[[i]]) &&
          !isTRUE(file.rename(backups[[i]], targets[[i]]))) {
        rollback <<- c(rollback, paste0("restore ", targets[[i]]))
      }
    }
    message <- conditionMessage(error)
    if (length(rollback)) {
      message <- paste0(message, "; rollback failed to: ",
                        paste(rollback, collapse = ", "))
    }
    .gsim_stop("Dataset triplet publication failed: ", message)
  })
  invisible(result)
}

.gsim_plink_dataset_create <- function(
  backend,
  metadata_backend,
  prefix,
  sample_metadata,
  overwrite = FALSE,
  buffer_variants = 64L,
  provenance = list()
) {
  if (!inherits(backend, "gsim_gbits_backend")) {
    .gsim_stop("backend must be created by .gsim_gbits_backend().")
  }
  validated_samples <- .gsim_plink_validate_samples(
    metadata_backend, sample_metadata
  )
  targets <- .gsim_plink_targets(prefix)
  overwrite <- .gsim_bed_sink_flag(overwrite, "overwrite")
  exists <- file.exists(targets)
  if (!overwrite && any(exists)) {
    .gsim_stop("PLINK destination exists and overwrite is disabled: ",
               paste(targets[exists], collapse = ", "))
  }
  if (overwrite && any(exists) && !all(exists)) {
    .gsim_stop("Overwrite requires an existing complete BED/BIM/FAM triplet.")
  }
  stage_prefix <- tempfile(
    pattern = ".gsim-stage-",
    tmpdir = dirname(targets[[1L]])
  )
  staged <- stats::setNames(paste0(stage_prefix, c(".bed", ".bim", ".fam")),
                            names(targets))
  state <- new.env(parent = emptyenv())
  state$bed <- .gsim_bed_sink_create(
    backend, staged[["bed"]], validated_samples$metadata$individual_id,
    overwrite = FALSE, buffer_variants = buffer_variants,
    provenance = provenance
  )
  state$backend <- backend
  state$metadata_backend <- metadata_backend
  state$sample_metadata <- validated_samples$metadata
  state$sample_pointer <- validated_samples$pointer
  state$targets <- targets
  state$staged <- staged
  state$overwrite <- overwrite
  state$variant_metadata <- list()
  state$provenance <- provenance
  state$finalized <- FALSE
  state$cancelled <- FALSE
  state$failed <- FALSE
  class(state) <- "gsim_plink_dataset"
  state
}

.gsim_plink_dataset_append <- function(dataset, chromosome, h1, h2,
                                       variant_metadata) {
  if (!inherits(dataset, "gsim_plink_dataset") ||
      !is.environment(dataset) || is.null(dataset$bed)) {
    .gsim_stop("dataset is not a valid experimental PLINK dataset sink.")
  }
  if (dataset$finalized || dataset$cancelled || dataset$failed) {
    .gsim_stop("Cannot append after dataset finalization, cancellation, or failure.")
  }
  metadata <- .gsim_plink_normalize_variants(chromosome, variant_metadata)
  packed_ids <- attr(h1, "variant_ids", exact = TRUE)
  if (!identical(metadata$variant_id, packed_ids) ||
      !identical(metadata$variant_id,
                 attr(h2, "variant_ids", exact = TRUE))) {
    .gsim_stop("Variant metadata order must exactly match packed H1/H2 IDs.")
  }
  # Validate the chromosome metadata before any BED record is appended.
  invisible(.gsim_gmat_variant_pointer(dataset$metadata_backend, metadata))
  .gsim_bed_sink_append(
    dataset$bed, chromosome, h1, h2, metadata$variant_id
  )
  dataset$variant_metadata[[length(dataset$variant_metadata) + 1L]] <- metadata
  invisible(dataset)
}

.gsim_plink_dataset_finalize <- function(
  dataset,
  .test_fail_stage = "none",
  .test_fail_publish_after = 0L
) {
  if (!inherits(dataset, "gsim_plink_dataset") || !is.environment(dataset)) {
    .gsim_stop("dataset is not a valid experimental PLINK dataset sink.")
  }
  if (dataset$finalized) .gsim_stop("PLINK dataset has already been finalized.")
  if (dataset$cancelled) .gsim_stop("A cancelled PLINK dataset cannot be finalized.")
  if (dataset$failed) .gsim_stop("A failed PLINK dataset cannot be finalized.")
  if (!length(dataset$variant_metadata)) {
    .gsim_stop("At least one chromosome must be appended before finalization.")
  }
  if (!identical(.test_fail_stage, "none") &&
      !identical(.test_fail_stage, "after_bed") &&
      !identical(.test_fail_stage, "after_bim")) {
    .gsim_stop(".test_fail_stage must be none, after_bed, or after_bim.")
  }
  success <- FALSE
  on.exit({
    if (!success) {
      dataset$failed <- TRUE
      if (!isTRUE(dataset$bed$finalized) && !isTRUE(dataset$bed$cancelled)) {
        try(.gsim_bed_sink_cancel(dataset$bed), silent = TRUE)
      }
      unlink(unname(dataset$staged), force = TRUE)
    }
  }, add = TRUE)

  variants <- do.call(rbind, dataset$variant_metadata)
  rownames(variants) <- NULL
  variant_pointer <- .gsim_gmat_variant_pointer(
    dataset$metadata_backend, variants
  )
  bed_manifest <- .gsim_bed_sink_finalize(dataset$bed)
  if (.test_fail_stage == "after_bed") {
    .gsim_stop("injected failure after BED completion")
  }
  bim_info <- .Call(
    C_gsim_gmat_write_bim, variant_pointer, enc2utf8(dataset$staged[["bim"]])
  )
  if (.test_fail_stage == "after_bim") {
    .gsim_stop("injected failure after BIM completion")
  }
  fam_info <- .Call(
    C_gsim_gmat_write_fam, dataset$sample_pointer,
    enc2utf8(dataset$staged[["fam"]])
  )

  expected_bed <- 3 + nrow(variants) *
    ceiling(nrow(dataset$sample_metadata) / 4)
  observed <- unname(file.info(dataset$staged)$size)
  if (anyNA(observed) || observed[[1L]] != expected_bed ||
      unname(bim_info[["record_count"]]) != nrow(variants) ||
      unname(fam_info[["record_count"]]) != nrow(dataset$sample_metadata)) {
    .gsim_stop("Staged PLINK triplet failed final count or byte validation.")
  }
  .gsim_dataset_publish(
    dataset$staged, dataset$targets, dataset$overwrite,
    .test_fail_after = .test_fail_publish_after
  )
  dataset$finalized <- TRUE
  success <- TRUE

  chromosome <- vapply(dataset$variant_metadata, function(x) x$chromosome[[1L]],
                       character(1))
  marker_count <- vapply(dataset$variant_metadata, nrow, integer(1))
  start <- cumsum(c(1L, head(marker_count, -1L)))
  manifest <- list(
    paths = dataset$targets,
    individual_count = nrow(dataset$sample_metadata),
    variant_count = nrow(variants),
    chromosomes = data.frame(
      chromosome = chromosome, marker_count = marker_count,
      start_variant = start, end_variant = start + marker_count - 1L,
      stringsAsFactors = FALSE
    ),
    expected_bed_bytes = expected_bed,
    observed_bed_bytes = unname(file.info(dataset$targets[["bed"]])$size),
    bim_record_count = unname(bim_info[["record_count"]]),
    fam_record_count = unname(fam_info[["record_count"]]),
    bim_bytes = unname(bim_info[["bytes_written"]]),
    fam_bytes = unname(fam_info[["bytes_written"]]),
    maximum_bim_record_bytes = unname(bim_info[["maximum_record_bytes"]]),
    maximum_fam_record_bytes = unname(fam_info[["maximum_record_bytes"]]),
    bed_conversion_buffer_bytes = bed_manifest$conversion_buffer_bytes,
    sample_ids = dataset$sample_metadata$individual_id,
    variant_ids = variants$variant_id,
    allele_orientation = "bit 1 = ALT = BIM A1; bit 0 = REF = BIM A2",
    backend = list(
      gbits = list(version = attr(dataset$backend, "gbits_version", exact = TRUE),
                   abi = attr(dataset$backend, "gbits_abi", exact = TRUE)),
      gmat = list(version = attr(dataset$metadata_backend, "gmat_version",
                                 exact = TRUE),
                  abi = attr(dataset$metadata_backend, "gmat_abi", exact = TRUE))
    ),
    provenance = dataset$provenance,
    publication_status = "published",
    transaction = "same-directory staging with backup-and-rollback publication"
  )
  class(manifest) <- c("gsim_plink_dataset_manifest", "list")
  manifest
}

.gsim_plink_dataset_cancel <- function(dataset) {
  if (!inherits(dataset, "gsim_plink_dataset") || !is.environment(dataset)) {
    .gsim_stop("dataset is not a valid experimental PLINK dataset sink.")
  }
  if (dataset$finalized) .gsim_stop("A finalized PLINK dataset cannot be cancelled.")
  if (dataset$cancelled) return(invisible(dataset))
  if (!isTRUE(dataset$bed$finalized) && !isTRUE(dataset$bed$cancelled)) {
    .gsim_bed_sink_cancel(dataset$bed)
  }
  unlink(c(unname(dataset$staged), paste0(unname(dataset$staged), ".backup")),
         force = TRUE)
  dataset$cancelled <- TRUE
  invisible(dataset)
}
