.ds_backend <- function() {
  gsim:::.gsim_gbits_backend(Sys.getenv("GSIM_GBITS_LIBRARY"))
}

.ds_metadata_backend <- function() {
  gsim:::.gsim_gmat_backend(Sys.getenv("GSIM_GMAT_LIBRARY"))
}

.ds_variants <- function(ids, chromosome = "chrZ", cm = seq_along(ids) - 1,
                         bp = seq_along(ids), alt = rep("G", length(ids)),
                         ref = rep("A", length(ids))) {
  data.frame(
    chromosome = rep(chromosome, length(ids)), variant_id = ids,
    genetic_position_cm = cm, base_pair_position = bp,
    alt = alt, ref = ref, bit1_allele = alt, bit0_allele = ref,
    stringsAsFactors = FALSE
  )
}

.ds_packed <- function(backend, h1, h2) {
  list(h1 = gsim:::.gsim_gbits_pack(backend, h1),
       h2 = gsim:::.gsim_gbits_pack(backend, h2))
}

.ds_owned_files <- function(directory) {
  list.files(directory, all.files = TRUE,
             pattern = "gsim-stage|gbits\\.tmp|gmat\\.tmp|\\.backup",
             full.names = TRUE)
}

.ds_founder_fixture <- function(backend, chromosome = "founder", markers = 13L,
                                n = 7L, seed = 1601L) {
  marker <- seq_len(markers)
  reference_h1 <- outer(seq_len(4L), marker,
                        function(i, j) (i + j + j %/% 2L) %% 2L)
  reference_h2 <- outer(seq_len(4L), marker,
                        function(i, j) (i * 3L + j * 5L + j %/% 3L) %% 2L)
  dimnames(reference_h1) <- dimnames(reference_h2) <- list(
    paste0("donor", seq_len(4L)), paste0(chromosome, "_v", marker)
  )
  args <- list(
    reference_haplotypes_h1 = reference_h1,
    reference_haplotypes_h2 = reference_h2,
    donor_population = c("A", "A", "B", "B"),
    ancestry_weights = c(A = 0.4, B = 0.6),
    N = c(A = 2, B = 2), Ne = c(A = 4, B = 6),
    rho = c(A = 0.8, B = 1.2),
    genetic_position = seq(0, 2, length.out = markers),
    mutation_age = rep(c(0.3, 3, 30, 1e8), length.out = markers),
    n = n, seed = seed, chromosome = rep(chromosome, markers)
  )
  list(
    raw = do.call(gsim:::.gsim_hapnest_founders,
                  c(args, list(return_genotypes = TRUE))),
    packed = do.call(gsim:::.gsim_hapnest_founders_packed_chromosome,
                     c(list(backend = backend), args,
                       list(return_genotypes = FALSE)))
  )
}

.ds_pedigree_fixture <- function(backend, chromosome = "pedigree",
                                 markers = 29L, seed = 812L) {
  tab <- data.frame(
    animal = c("S", "D", "D2", "S2", "C", "FS", "PHS", "MHS", "G", "L"),
    sire = c(NA, NA, NA, NA, "S", "S", "S", "S2", "C", NA),
    dam = c(NA, NA, NA, NA, "D", "D", "D2", "D", "D2", NA),
    stringsAsFactors = FALSE
  )
  pedigree <- structure(list(
    pedigree = tab[c(10, 6, 3, 9, 1, 7, 4, 5, 2, 8), , drop = FALSE],
    canonical_order = tab$animal,
    external_order = tab$animal[c(10, 6, 3, 9, 1, 7, 4, 5, 2, 8)]
  ), class = "gsim_pedigree")
  founders <- c("S", "D", "D2", "S2", "L")
  variants <- paste0(chromosome, "_v", seq_len(markers))
  h1 <- outer(seq_along(founders), seq_len(markers),
              function(i, j) (i + j + j %/% 4L) %% 2L)
  h2 <- outer(seq_along(founders), seq_len(markers),
              function(i, j) (i * 3L + j * 5L + j %/% 3L + 1L) %% 2L)
  dimnames(h1) <- dimnames(h2) <- list(founders, variants)
  args <- list(
    pedigree = pedigree, founder_haplotypes = list(h1 = h1, h2 = h2),
    chromosome = rep(chromosome, markers),
    genetic_position = seq(0, 10, length.out = markers), seed = seed,
    return_crossovers = TRUE
  )
  list(
    pedigree = pedigree,
    raw = do.call(gsim:::.gsim_pedigree_genotypes,
                  c(args, list(return_genotypes = TRUE))),
    packed = do.call(gsim:::.gsim_pedigree_genotypes_packed_chromosome,
                     c(list(backend = backend), args,
                       list(return_genotypes = FALSE)))
  )
}

testthat::test_that("dataset freezes exact BED, BIM, FAM, and allele bytes", {
  backend <- .ds_backend()
  metadata_backend <- .ds_metadata_backend()
  root <- tempfile("plink-exact-")
  dir.create(root)
  directory <- file.path(root, "directory with spaces")
  dir.create(directory)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ids <- c("S", "D", "C", "L")
  samples <- gsim:::.gsim_plink_sample_metadata(
    ids, family_id = rep("F", 4L), paternal_id = c("0", "0", "S", "0"),
    maternal_id = c("0", "0", "D", "0"), sex = c(1L, 2L, 0L, 0L)
  )
  h1 <- matrix(as.raw(c(0, 0, 1, 1)), 4L, 1L,
               dimnames = list(ids, "v1"))
  h2 <- matrix(as.raw(c(0, 1, 0, 1)), 4L, 1L, dimnames = dimnames(h1))
  packed <- .ds_packed(backend, h1, h2)
  dataset <- gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, "dæta set"), samples,
    provenance = list(seed = 1L, source = "hand fixture")
  )
  gsim:::.gsim_plink_dataset_append(
    dataset, "chrZ", packed$h1, packed$h2,
    .ds_variants("v1", "chrZ", 0, 101, "G", "A")
  )
  manifest <- gsim:::.gsim_plink_dataset_finalize(dataset)
  testthat::expect_identical(
    readBin(manifest$paths[["bed"]], "raw", n = 99L),
    as.raw(c(0x6c, 0x1b, 0x01, 0x2b))
  )
  testthat::expect_identical(
    readBin(manifest$paths[["bim"]], "raw", n = 99L),
    charToRaw("chrZ\tv1\t0\t101\tG\tA\n")
  )
  testthat::expect_identical(
    readBin(manifest$paths[["fam"]], "raw", n = 999L),
    charToRaw(paste0(
      "F\tS\t0\t0\t1\t-9\n", "F\tD\t0\t0\t2\t-9\n",
      "F\tC\tS\tD\t0\t-9\n", "F\tL\t0\t0\t0\t-9\n"
    ))
  )
  decoded <- gsim:::.gsim_gbits_bed_read_all(
    backend, manifest$paths[["bed"]], 4L, 1L, ids, "v1"
  )
  testthat::expect_identical(decoded, matrix(c(0L, 1L, 1L, 2L), 4L, 1L,
                                             dimnames = list(ids, "v1")))
  testthat::expect_identical(manifest$allele_orientation,
    "bit 1 = ALT = BIM A1; bit 0 = REF = BIM A2")
  testthat::expect_identical(manifest$backend$gbits,
                             list(version = "0.20.0", abi = 4L))
  testthat::expect_identical(manifest$backend$gmat,
                             list(version = "0.4.0", abi = 0L))
  testthat::expect_identical(manifest$publication_status, "published")
})

testthat::test_that("BIM validation rejects ambiguous identity, maps, and alleles", {
  backend <- .ds_backend()
  metadata_backend <- .ds_metadata_backend()
  directory <- tempfile("plink-bim-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  ids <- c("a", "b")
  samples <- gsim:::.gsim_plink_sample_metadata(ids)
  h <- matrix(as.raw(c(0, 1, 1, 0)), 2L, 2L,
              dimnames = list(ids, c("v1", "v2")))
  packed <- .ds_packed(backend, h, h)
  make <- function(name) gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, name), samples
  )
  invalid <- list(
    duplicate = .ds_variants(c("v1", "v1")),
    empty = .ds_variants(c("v1", "")),
    missing = .ds_variants(c("v1", NA_character_)),
    bad_cm = .ds_variants(c("v1", "v2"), cm = c(0, -1)),
    unsorted_cm = .ds_variants(c("v1", "v2"), cm = c(1, 0)),
    bad_bp = .ds_variants(c("v1", "v2"), bp = c(1, 0)),
    unsorted_bp = .ds_variants(c("v1", "v2"), bp = c(2, 1)),
    identical_alleles = .ds_variants(c("v1", "v2"),
                                     alt = c("G", "A"), ref = c("A", "A")),
    multibase = .ds_variants(c("v1", "v2"), alt = c("G", "DEL")),
    symbolic = .ds_variants(c("v1", "v2"), alt = c("G", "<DEL>"))
  )
  for (name in names(invalid)) {
    dataset <- make(name)
    testthat::expect_error(gsim:::.gsim_plink_dataset_append(
      dataset, "chrZ", packed$h1, packed$h2, invalid[[name]]
    ))
    gsim:::.gsim_plink_dataset_cancel(dataset)
  }
  reversed <- .ds_variants(c("v1", "v2"))
  reversed$bit1_allele <- reversed$ref
  reversed$bit0_allele <- reversed$alt
  dataset <- make("reversed")
  testthat::expect_error(gsim:::.gsim_plink_dataset_append(
    dataset, "chrZ", packed$h1, packed$h2, reversed
  ), "allele reversal is not implicit")
  gsim:::.gsim_plink_dataset_cancel(dataset)

  mismatched <- .ds_variants(c("v2", "v1"))
  dataset <- make("order")
  testthat::expect_error(gsim:::.gsim_plink_dataset_append(
    dataset, "chrZ", packed$h1, packed$h2, mismatched
  ), "order")
  gsim:::.gsim_plink_dataset_cancel(dataset)
})

testthat::test_that("FAM validation preserves pedigree roles and rejects ambiguity", {
  backend <- .ds_backend()
  metadata_backend <- .ds_metadata_backend()
  directory <- tempfile("plink-fam-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  valid <- gsim:::.gsim_plink_sample_metadata(
    c("S", "D", "D2", "C", "FS", "PHS", "G", "L"),
    paternal_id = c("0", "0", "0", "S", "S", "S", "C", "0"),
    maternal_id = c("0", "0", "0", "D", "D", "D2", "D2", "0"),
    sex = c(1L, 2L, 2L, 1L, 0L, 0L, 0L, 0L)
  )
  make <- function(metadata, name) gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, name), metadata
  )
  testthat::expect_s3_class(make(valid, "valid"), "gsim_plink_dataset")

  cases <- list(
    duplicate = transform(valid, individual_id = replace(individual_id, 2L, "S")),
    missing = transform(valid, individual_id = replace(individual_id, 2L, NA)),
    self = transform(valid, paternal_id = replace(paternal_id, 4L, "C")),
    partial = transform(valid, maternal_id = replace(maternal_id, 4L, "0")),
    absent = transform(valid, paternal_id = replace(paternal_id, 4L, "absent")),
    parent_after = valid[c(4L, 1L, 2L, 3L, 5L, 6L, 7L, 8L), ],
    bad_sex = transform(valid, sex = replace(sex, 1L, 2L)),
    phenotype = transform(valid, phenotype = replace(phenotype, 1L, 1.5))
  )
  for (name in names(cases)) {
    testthat::expect_error(make(cases[[name]], name))
  }

  ids <- valid$individual_id
  h <- matrix(as.raw(rep(c(0, 1), length.out = length(ids))), length(ids), 1L,
              dimnames = list(ids, "v"))
  permuted <- valid[c(2L, 1L, 3:8), , drop = FALSE]
  dataset <- make(permuted, "permuted")
  testthat::expect_error(gsim:::.gsim_plink_dataset_append(
    dataset, "x", gsim:::.gsim_gbits_pack(backend, h),
    gsim:::.gsim_gbits_pack(backend, h), .ds_variants("v", "x")
  ), "sample order")
  gsim:::.gsim_plink_dataset_cancel(dataset)
})

testthat::test_that("multi-chromosome BIM and BED order are literal and exact", {
  backend <- .ds_backend()
  metadata_backend <- .ds_metadata_backend()
  directory <- tempfile("plink-order-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  ids <- paste0("i", 1:5)
  dataset <- gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, "ordered"),
    gsim:::.gsim_plink_sample_metadata(ids), buffer_variants = 1L
  )
  labels <- c("chrZ", "01", "CHR1")
  marker_counts <- c(2L, 1L, 3L)
  all_genotypes <- list()
  expected_bim <- character()
  for (block in seq_along(labels)) {
    variants <- paste0(labels[[block]], "_v", seq_len(marker_counts[[block]]))
    h1 <- outer(seq_along(ids), seq_along(variants),
                function(i, j) (i + j + block) %% 2L)
    h2 <- outer(seq_along(ids), seq_along(variants),
                function(i, j) (i * 3L + j + block) %% 2L)
    dimnames(h1) <- dimnames(h2) <- list(ids, variants)
    metadata <- .ds_variants(
      variants, labels[[block]], c(0, seq_len(length(variants) - 1L) / 8),
      block * 100L + seq_along(variants),
      rep(c("G", "T"), length.out = length(variants)),
      rep(c("A", "C"), length.out = length(variants))
    )
    packed <- .ds_packed(backend, h1, h2)
    gsim:::.gsim_plink_dataset_append(
      dataset, labels[[block]], packed$h1, packed$h2, metadata
    )
    all_genotypes[[block]] <- h1 + h2
    cm_text <- vapply(metadata$genetic_position_cm, function(value) {
      text <- formatC(value, format = "f", digits = 15)
      text <- sub("0+$", "", text)
      sub("\\.$", "", text)
    }, character(1))
    expected_bim <- c(expected_bim, paste(
      metadata$chromosome, metadata$variant_id,
      cm_text,
      metadata$base_pair_position, metadata$alt, metadata$ref, sep = "\t"
    ))
    rm(packed, h1, h2)
    gc(FALSE)
  }
  manifest <- gsim:::.gsim_plink_dataset_finalize(dataset)
  decoded <- gsim:::.gsim_gbits_bed_read_all(
    backend, manifest$paths[["bed"]], length(ids), sum(marker_counts), ids,
    manifest$variant_ids
  )
  expected <- do.call(cbind, all_genotypes)
  testthat::expect_identical(decoded,
    matrix(as.integer(expected), nrow(expected), dimnames = dimnames(expected)))
  testthat::expect_identical(readLines(manifest$paths[["bim"]]), expected_bim)
  testthat::expect_identical(manifest$chromosomes$chromosome, labels)
  testthat::expect_identical(manifest$chromosomes$marker_count, marker_counts)
  testthat::expect_identical(manifest$chromosomes$start_variant, c(1L, 3L, 4L))
  expected_size <- 3 + sum(marker_counts) * ceiling(length(ids) / 4)
  testthat::expect_identical(manifest$expected_bed_bytes, expected_size)
  testthat::expect_identical(manifest$observed_bed_bytes, expected_size)
})

testthat::test_that("triplet transaction stages, cancels, overwrites, and rolls back", {
  backend <- .ds_backend()
  metadata_backend <- .ds_metadata_backend()
  root <- tempfile("plink-transaction-")
  dir.create(root)
  directory <- file.path(root, "directory with spaces")
  dir.create(directory)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  sentinel <- file.path(directory, "unrelated.keep")
  writeBin(charToRaw("keep"), sentinel)
  ids <- c("a", "b")
  samples <- gsim:::.gsim_plink_sample_metadata(ids)
  h1 <- matrix(as.raw(c(0, 1)), 2L, 1L, dimnames = list(ids, "v"))
  h2 <- matrix(as.raw(c(1, 1)), 2L, 1L, dimnames = dimnames(h1))
  packed <- .ds_packed(backend, h1, h2)
  metadata <- .ds_variants("v", "u-β", 0.25, 9, "T", "C")
  prefix <- file.path(directory, "data set æ")
  targets <- paste0(prefix, c(".bed", ".bim", ".fam"))
  make <- function(overwrite = FALSE) {
    dataset <- gsim:::.gsim_plink_dataset_create(
      backend, metadata_backend, prefix, samples, overwrite = overwrite
    )
    gsim:::.gsim_plink_dataset_append(
      dataset, "u-β", packed$h1, packed$h2, metadata
    )
    dataset
  }

  cancelled <- make()
  gsim:::.gsim_plink_dataset_cancel(cancelled)
  testthat::expect_false(any(file.exists(targets)))
  testthat::expect_length(.ds_owned_files(directory), 0L)
  testthat::expect_identical(readBin(sentinel, "raw", n = 9L),
                             charToRaw("keep"))

  before_bed <- gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, "invalid"), samples
  )
  bad <- metadata
  bad$variant_id <- "different"
  testthat::expect_error(gsim:::.gsim_plink_dataset_append(
    before_bed, "u-β", packed$h1, packed$h2, bad
  ), "order")
  gsim:::.gsim_plink_dataset_cancel(before_bed)

  after_bed <- make()
  testthat::expect_error(gsim:::.gsim_plink_dataset_finalize(
    after_bed, .test_fail_stage = "after_bed"
  ), "after BED")
  testthat::expect_false(any(file.exists(targets)))
  testthat::expect_length(.ds_owned_files(directory), 0L)

  after_bim <- make()
  testthat::expect_error(gsim:::.gsim_plink_dataset_finalize(
    after_bim, .test_fail_stage = "after_bim"
  ), "after BIM")
  testthat::expect_false(any(file.exists(targets)))
  testthat::expect_length(.ds_owned_files(directory), 0L)

  manifest <- gsim:::.gsim_plink_dataset_finalize(make())
  first <- lapply(targets, readBin, what = "raw", n = 999L)
  testthat::expect_error(make(), "overwrite is disabled")

  replacement <- make(overwrite = TRUE)
  testthat::expect_identical(lapply(targets, readBin, what = "raw", n = 999L),
                             first)
  second <- gsim:::.gsim_plink_dataset_finalize(replacement)
  testthat::expect_identical(second$publication_status, "published")

  old <- list(charToRaw("old-bed"), charToRaw("old-bim"), charToRaw("old-fam"))
  Map(writeBin, old, targets)
  failing <- make(overwrite = TRUE)
  testthat::expect_error(gsim:::.gsim_plink_dataset_finalize(
    failing, .test_fail_publish_after = 1L
  ), "publication failed")
  testthat::expect_identical(lapply(targets, readBin, what = "raw", n = 999L),
                             old)
  testthat::expect_length(.ds_owned_files(directory), 0L)
  testthat::expect_false(isTRUE(failing$finalized))
  testthat::expect_true(isTRUE(failing$failed))
  testthat::expect_identical(readBin(sentinel, "raw", n = 9L),
                             charToRaw("keep"))

  unlink(targets[[3L]])
  testthat::expect_error(make(overwrite = TRUE), "complete BED/BIM/FAM")
  testthat::expect_true(file.exists(targets[[1L]]) && file.exists(targets[[2L]]))
  testthat::expect_false(file.exists(targets[[3L]]))
  invisible(manifest)
})

testthat::test_that("founder and multigeneration datasets decode to both oracles", {
  backend <- .ds_backend()
  metadata_backend <- .ds_metadata_backend()
  directory <- tempfile("plink-parity-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)

  founder <- .ds_founder_fixture(backend)
  founder_data <- gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, "founder"),
    gsim:::.gsim_plink_sample_metadata(founder$packed$sample_ids)
  )
  gsim:::.gsim_plink_dataset_append(
    founder_data, founder$packed$chromosome, founder$packed$h1,
    founder$packed$h2,
    .ds_variants(founder$packed$variant_ids, founder$packed$chromosome,
                 seq(0, 200, length.out = length(founder$packed$variant_ids)),
                 seq_along(founder$packed$variant_ids))
  )
  founder_manifest <- gsim:::.gsim_plink_dataset_finalize(founder_data)
  founder_decoded <- gsim:::.gsim_gbits_bed_read_all(
    backend, founder_manifest$paths[["bed"]],
    length(founder$packed$sample_ids), length(founder$packed$variant_ids),
    founder_manifest$sample_ids, founder_manifest$variant_ids
  )
  testthat::expect_identical(founder_decoded,
    matrix(as.integer(founder$raw$genotypes), nrow(founder$raw$genotypes),
           dimnames = dimnames(founder$raw$genotypes)))
  testthat::expect_identical(founder_decoded,
    matrix(as.integer(gsim:::.gsim_gbits_decode_genotypes(
      founder$packed$h1, founder$packed$h2)), nrow(founder_decoded),
      dimnames = dimnames(founder_decoded)))

  pedigree <- .ds_pedigree_fixture(backend)
  sample_metadata <- gsim:::.gsim_plink_pedigree_metadata(
    pedigree$pedigree,
    sex = c(1L, 2L, 2L, 1L, 1L, 0L, 0L, 0L, 0L, 0L)
  )
  pedigree_data <- gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, "pedigree"), sample_metadata
  )
  pedigree_chromosome <- pedigree$packed$chromosome_blocks$chromosome[[1L]]
  gsim:::.gsim_plink_dataset_append(
    pedigree_data, pedigree_chromosome, pedigree$packed$h1,
    pedigree$packed$h2,
    .ds_variants(pedigree$packed$variant_ids, pedigree_chromosome,
                 seq(0, 1000, length.out = length(pedigree$packed$variant_ids)),
                 seq_along(pedigree$packed$variant_ids))
  )
  pedigree_manifest <- gsim:::.gsim_plink_dataset_finalize(pedigree_data)
  decoded <- gsim:::.gsim_gbits_bed_read_all(
    backend, pedigree_manifest$paths[["bed"]],
    length(pedigree$packed$sample_ids), length(pedigree$packed$variant_ids),
    pedigree_manifest$sample_ids, pedigree_manifest$variant_ids
  )
  expected <- matrix(as.integer(pedigree$raw$genotypes),
                     nrow(pedigree$raw$genotypes),
                     dimnames = dimnames(pedigree$raw$genotypes))
  testthat::expect_identical(decoded, expected)
  testthat::expect_identical(decoded,
    matrix(as.integer(gsim:::.gsim_gbits_decode_genotypes(
      pedigree$packed$h1, pedigree$packed$h2)), nrow(decoded),
      dimnames = dimnames(decoded)))
  testthat::expect_identical(sum(decoded == -9L), 0L)
  fam <- read.delim(pedigree_manifest$paths[["fam"]], header = FALSE,
                    colClasses = "character")
  testthat::expect_identical(fam[[2L]], pedigree$packed$sample_ids)
  testthat::expect_identical(fam[[3L]], sample_metadata$paternal_id)
  testthat::expect_identical(fam[[4L]], sample_metadata$maternal_id)

  tab <- pedigree$pedigree$pedigree
  tab <- tab[match(pedigree$raw$sample_ids, tab$animal), , drop = FALSE]
  inconsistency <- 0L
  for (i in which(!is.na(tab$sire))) {
    sire <- match(tab$sire[[i]], tab$animal)
    dam <- match(tab$dam[[i]], tab$animal)
    # Marker-wise maternal and paternal checks, deliberately independent of
    # the production gamete materializer.
    inconsistency <- inconsistency + sum(vapply(seq_len(ncol(expected)), function(j) {
      !(pedigree$raw$h1[i, j] %in%
          c(pedigree$raw$h1[sire, j], pedigree$raw$h2[sire, j])) ||
        !(pedigree$raw$h2[i, j] %in%
            c(pedigree$raw$h1[dam, j], pedigree$raw$h2[dam, j]))
    }, logical(1)))
  }
  testthat::expect_identical(inconsistency, 0L)
  testthat::expect_true(all(c("FS", "PHS", "MHS", "G", "L") %in%
                              pedigree_manifest$sample_ids))
})

testthat::test_that("dataset memory accounting remains chromosome-wise", {
  backend <- .ds_backend()
  metadata_backend <- .ds_metadata_backend()
  directory <- tempfile("plink-memory-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  individuals <- 257L
  markers <- c(127L, 65L, 1L)
  ids <- paste0("i", seq_len(individuals))
  dataset <- gsim:::.gsim_plink_dataset_create(
    backend, metadata_backend, file.path(directory, "bounded"),
    gsim:::.gsim_plink_sample_metadata(ids), buffer_variants = 7L
  )
  peak_packed <- 0
  elapsed <- system.time(for (block in seq_along(markers)) {
    variants <- paste0("c", block, "v", seq_len(markers[[block]]))
    h1 <- matrix(as.raw(rep(c(0, 1), length.out = individuals * markers[[block]])),
                 individuals, markers[[block]], dimnames = list(ids, variants))
    h2 <- matrix(as.raw(rep(c(1, 0, 0), length.out = individuals * markers[[block]])),
                 individuals, markers[[block]], dimnames = list(ids, variants))
    packed <- .ds_packed(backend, h1, h2)
    packed_bytes <- unname(gsim:::.gsim_gbits_info(packed$h1)[[4L]] +
                               gsim:::.gsim_gbits_info(packed$h2)[[4L]])
    peak_packed <- max(peak_packed, packed_bytes)
    gsim:::.gsim_plink_dataset_append(
      dataset, paste0("c", block), packed$h1, packed$h2,
      .ds_variants(variants, paste0("c", block),
                   seq(0, 1, length.out = markers[[block]]),
                   seq_len(markers[[block]]))
    )
    rm(h1, h2, packed)
    gc(FALSE)
  })[["elapsed"]]
  manifest <- gsim:::.gsim_plink_dataset_finalize(dataset)
  expected_peak <- 2 * 8 * max(markers) * ceiling(individuals / 64)
  expected_bed <- 3 + sum(markers) * ceiling(individuals / 4)
  expected_buffer <- 7 * ceiling(individuals / 4)
  testthat::expect_identical(peak_packed, expected_peak)
  testthat::expect_identical(manifest$observed_bed_bytes, expected_bed)
  testthat::expect_identical(manifest$bed_conversion_buffer_bytes,
                             expected_buffer)
  testthat::expect_lte(manifest$maximum_bim_record_bytes, 128)
  testthat::expect_lte(manifest$maximum_fam_record_bytes, 128)
  testthat::expect_false("genotypes" %in% names(manifest))
  testthat::expect_true(is.finite(elapsed) && elapsed < 30)
})
