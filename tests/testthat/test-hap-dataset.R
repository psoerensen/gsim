.hap_backend <- function() {
  gsim:::.gsim_packed_backend()
}

.hap_metadata_backend <- function() {
  gsim:::.gsim_metadata_backend()
}

.hap_variants <- function(chromosome, ids, offset = 0L) {
  n <- length(ids)
  data.frame(
    chromosome = rep(chromosome, n), variant_id = ids,
    genetic_position_cm = seq(0, by = 0.125, length.out = n),
    base_pair_position = offset + seq_len(n),
    alt = rep(c("C", "G"), length.out = n),
    ref = rep(c("A", "T"), length.out = n),
    bit1_allele = rep(c("C", "G"), length.out = n),
    bit0_allele = rep(c("A", "T"), length.out = n),
    stringsAsFactors = FALSE
  )
}

.hap_pack_pair <- function(backend, h1, h2) {
  list(h1 = gsim:::.gsim_packed_pack(backend, h1),
       h2 = gsim:::.gsim_packed_pack(backend, h2))
}

.hap_mendelian_errors <- function(h1, h2, alignment) {
  errors <- 0L
  for (i in which(!alignment$founder)) {
    sire <- match(alignment$sire[[i]], alignment$animal)
    dam <- match(alignment$dam[[i]], alignment$animal)
    errors <- errors + sum(!(h1[i, ] == h1[sire, ] | h1[i, ] == h2[sire, ]))
    errors <- errors + sum(!(h2[i, ] == h1[dam, ] | h2[i, ] == h2[dam, ]))
  }
  errors
}

testthat::test_that("HAP v1 writes exact header, ranges, and packed round trips", {
  backend <- .hap_backend()
  metadata_backend <- .hap_metadata_backend()
  root <- tempfile("hap-v1-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  samples <- paste0("s", seq_len(65L))
  sample_metadata <- gsim:::.gsim_plink_sample_metadata(samples)
  set.seed(2468)
  rng_before <- .Random.seed
  chromosomes <- list(`chr:10` = 3L, `01` = 1L, `Z-alt` = 5L)
  pairs <- list()
  dataset <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "dÃ¦ta set"), sample_metadata,
    provenance = list(seed = 17L, source = "hand packed fixture")
  )
  offset <- 0L
  for (label in names(chromosomes)) {
    m <- chromosomes[[label]]
    ids <- paste0(label, "_v", seq_len(m))
    h1 <- outer(seq_len(65L), seq_len(m), function(i, j) (i + j) %% 2L)
    h2 <- outer(seq_len(65L), seq_len(m), function(i, j) (i * 3L + j + 1L) %% 2L)
    dimnames(h1) <- dimnames(h2) <- list(samples, ids)
    pairs[[label]] <- .hap_pack_pair(backend, h1, h2)
    gsim:::.gsim_hap_dataset_append(
      dataset, label, pairs[[label]]$h1, pairs[[label]]$h2,
      .hap_variants(label, ids, offset)
    )
    offset <- offset + m
  }
  manifest <- gsim:::.gsim_hap_dataset_finalize(dataset)
  testthat::expect_identical(
    readBin(manifest$paths[["hap"]], "raw", n = 4L),
    as.raw(c(0x48, 0x41, 0x50, 0x01))
  )
  expected_payload <- 2 * sum(unlist(chromosomes)) * ceiling(65 / 64) * 8
  expected_size <- 64 + expected_payload + 48 * length(chromosomes)
  testthat::expect_equal(manifest$expected_hap_bytes, expected_size)
  testthat::expect_equal(manifest$observed_hap_bytes, expected_size)
  testthat::expect_equal(manifest$header_bytes, 64)
  testthat::expect_equal(manifest$chromosome_table_bytes, 144)
  testthat::expect_equal(manifest$packed_data_bytes, expected_payload)
  testthat::expect_identical(manifest$implementation$engine,
                             "gsim private native backend")
  testthat::expect_match(manifest$implementation$packed_origin, "089bf1e")
  testthat::expect_match(manifest$implementation$metadata_origin, "33d6751")

  reader <- gsim:::.gsim_hap_dataset_open(
    backend, metadata_backend, file.path(root, "dÃ¦ta set"))
  inspected <- gsim:::.gsim_hap_dataset_inspect(reader)
  testthat::expect_identical(inspected$chromosomes$chromosome,
                             names(chromosomes))
  for (label in rev(names(chromosomes))) {
    loaded <- gsim:::.gsim_hap_dataset_load_chromosome(reader, label)
    testthat::expect_identical(
      gsim:::.gsim_packed_unpack(loaded$h1),
      gsim:::.gsim_packed_unpack(pairs[[label]]$h1))
    testthat::expect_identical(
      gsim:::.gsim_packed_unpack(loaded$h2),
      gsim:::.gsim_packed_unpack(pairs[[label]]$h2))
    info_loaded <- gsim:::.gsim_packed_info(loaded$h1)
    info_original <- gsim:::.gsim_packed_info(pairs[[label]]$h1)
    for (marker in seq_len(chromosomes[[label]])) {
      for (word in seq_len(info_loaded[[3L]])) {
        testthat::expect_identical(
          gsim:::.gsim_packed_word(loaded$h1, marker, word),
          gsim:::.gsim_packed_word(pairs[[label]]$h1, marker, word))
      }
    }
    testthat::expect_identical(info_loaded, info_original)
  }
  retained <- gsim:::.gsim_hap_dataset_load_chromosome(reader, "01")
  gsim:::.gsim_hap_dataset_close(reader)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(retained$h1),
                             gsim:::.gsim_packed_unpack(pairs[["01"]]$h1))
  testthat::expect_identical(.Random.seed, rng_before)
})

testthat::test_that("HAP metadata alignment and lifecycle failures are strict", {
  backend <- .hap_backend()
  metadata_backend <- .hap_metadata_backend()
  root <- tempfile("hap-invalid-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  ids <- c("a", "b", "c")
  variants <- c("v1", "v2")
  h <- matrix(as.raw(c(0, 1, 0, 1, 0, 1)), 3L, 2L,
              dimnames = list(ids, variants))
  phases <- .hap_pack_pair(backend, h, h)
  samples <- gsim:::.gsim_plink_sample_metadata(ids)

  cancelled <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "cancelled"), samples)
  gsim:::.gsim_hap_dataset_cancel(cancelled)
  testthat::expect_false(any(file.exists(
    paste0(file.path(root, "cancelled"), c(".hap", ".bim", ".fam")))))

  failed <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "failed"), samples)
  gsim:::.gsim_hap_dataset_append(
    failed, "x", phases$h1, phases$h2, .hap_variants("x", variants))
  testthat::expect_error(
    gsim:::.gsim_hap_dataset_finalize(failed, .test_fail_stage = "after_hap"),
    "injected failure")
  testthat::expect_false(any(file.exists(
    paste0(file.path(root, "failed"), c(".hap", ".bim", ".fam")))))

  good <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "good"), samples)
  gsim:::.gsim_hap_dataset_append(
    good, "x", phases$h1, phases$h2, .hap_variants("x", variants))
  manifest <- gsim:::.gsim_hap_dataset_finalize(good)
  testthat::expect_error(gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "good"), samples), "exists")

  old_bytes <- lapply(manifest$paths, readBin, what = "raw", n = 10000L)
  replacement <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "good"), samples,
    overwrite = TRUE)
  gsim:::.gsim_hap_dataset_append(
    replacement, "x", phases$h1, phases$h2, .hap_variants("x", variants))
  testthat::expect_error(gsim:::.gsim_hap_dataset_finalize(
    replacement, .test_fail_publish_after = 1L), "publication failed")
  testthat::expect_true(all(file.exists(manifest$paths)))
  testthat::expect_identical(
    lapply(manifest$paths, readBin, what = "raw", n = 10000L), old_bytes)

  replacement <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "good"), samples,
    overwrite = TRUE)
  gsim:::.gsim_hap_dataset_append(
    replacement, "x", phases$h1, phases$h2, .hap_variants("x", variants))
  replaced <- gsim:::.gsim_hap_dataset_finalize(replacement)
  testthat::expect_identical(replaced$publication_status, "published")

  # Count mismatch is detectable and rejected. Same-shape metadata substitution
  # is intentionally not claimed detectable because HAP v1 has no checksums.
  fam_bytes <- readBin(replaced$paths[["fam"]], "raw", n = 10000L)
  writeBin(c(fam_bytes, charToRaw("F\textra\t0\t0\t0\t-9\n")),
           replaced$paths[["fam"]])
  testthat::expect_error(gsim:::.gsim_hap_dataset_open(
    backend, metadata_backend, file.path(root, "good")), "align")
})

testthat::test_that("founder and pedigree packed phases survive HAP and emit exact BED dosages", {
  backend <- .hap_backend()
  metadata_backend <- .hap_metadata_backend()
  root <- tempfile("hap-simulation-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  markers <- 13L
  reference_h1 <- outer(1:4, 1:markers, function(i, j) (i + j) %% 2L)
  reference_h2 <- outer(1:4, 1:markers, function(i, j) (i * 3L + j + 1L) %% 2L)
  dimnames(reference_h1) <- dimnames(reference_h2) <- list(
    paste0("d", 1:4), paste0("f", 1:markers))
  founder_args <- list(
    reference_haplotypes_h1 = reference_h1,
    reference_haplotypes_h2 = reference_h2,
    donor_population = c("A", "A", "B", "B"),
    ancestry_weights = c(A = 0.4, B = 0.6), N = c(A = 2, B = 2),
    Ne = c(A = 4, B = 6), rho = c(A = 0.8, B = 1.2),
    genetic_position = seq(0, 2, length.out = markers),
    mutation_age = rep(c(0.2, 2, 20, 1e8), length.out = markers),
    n = 6L, seed = 101L, chromosome = rep("founder", markers))
  raw_founder <- do.call(gsim:::.gsim_hapnest_founders,
                         c(founder_args, list(return_genotypes = TRUE)))
  packed_founder <- do.call(gsim:::.gsim_hapnest_founders_packed_chromosome,
     c(list(backend = backend), founder_args, list(return_genotypes = FALSE)))
  founder_samples <- gsim:::.gsim_plink_sample_metadata(packed_founder$sample_ids)
  founder_dataset <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "founders"), founder_samples)
  gsim:::.gsim_hap_dataset_append(
    founder_dataset, "founder", packed_founder$h1, packed_founder$h2,
    .hap_variants("founder", packed_founder$variant_ids))
  gsim:::.gsim_hap_dataset_finalize(founder_dataset)
  founder_reader <- gsim:::.gsim_hap_dataset_open(
    backend, metadata_backend, file.path(root, "founders"))
  reloaded_founder <- gsim:::.gsim_hap_dataset_load_chromosome(
    founder_reader, "founder")
  testthat::expect_identical(gsim:::.gsim_packed_unpack(reloaded_founder$h1),
                             raw_founder$h1)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(reloaded_founder$h2),
                             raw_founder$h2)
  gsim:::.gsim_hap_dataset_close(founder_reader)

  tab <- data.frame(
    animal = c("S", "D", "D2", "S2", "C", "FS", "PHS", "MHS", "G", "L"),
    sire = c(NA, NA, NA, NA, "S", "S", "S", "S2", "C", NA),
    dam = c(NA, NA, NA, NA, "D", "D", "D2", "D", "D2", NA),
    stringsAsFactors = FALSE)
  pedigree <- structure(list(pedigree = tab[c(10,6,3,9,1,7,4,5,2,8), ],
                              canonical_order = tab$animal,
                              external_order = tab$animal[c(10,6,3,9,1,7,4,5,2,8)]),
                         class = "gsim_pedigree")
  founder_ids <- c("S", "D", "D2", "S2", "L")
  pvars <- paste0("p", seq_len(25L))
  fh1 <- outer(seq_along(founder_ids), seq_along(pvars),
               function(i, j) (i + j + j %/% 4L) %% 2L)
  fh2 <- outer(seq_along(founder_ids), seq_along(pvars),
               function(i, j) (i * 3L + j * 5L + 1L) %% 2L)
  dimnames(fh1) <- dimnames(fh2) <- list(founder_ids, pvars)
  meiosis_args <- list(
    pedigree = pedigree, founder_haplotypes = list(h1 = fh1, h2 = fh2),
    chromosome = rep("ped", length(pvars)),
    genetic_position = seq(0, 8, length.out = length(pvars)), seed = 909L,
    return_crossovers = TRUE)
  raw_pedigree <- do.call(gsim:::.gsim_pedigree_genotypes,
                          c(meiosis_args, list(return_genotypes = TRUE)))
  packed_pedigree <- do.call(gsim:::.gsim_pedigree_genotypes_packed_chromosome,
     c(list(backend = backend), meiosis_args, list(return_genotypes = FALSE)))
  pedigree_samples <- gsim:::.gsim_plink_pedigree_metadata(pedigree)
  pedigree_dataset <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, file.path(root, "pedigree"), pedigree_samples)
  gsim:::.gsim_hap_dataset_append(
    pedigree_dataset, "ped", packed_pedigree$h1, packed_pedigree$h2,
    .hap_variants("ped", pvars))
  gsim:::.gsim_hap_dataset_finalize(pedigree_dataset)
  pedigree_reader <- gsim:::.gsim_hap_dataset_open(
    backend, metadata_backend, file.path(root, "pedigree"))
  reloaded <- gsim:::.gsim_hap_dataset_load_chromosome(pedigree_reader, "ped")
  loaded_h1 <- gsim:::.gsim_packed_unpack(reloaded$h1)
  loaded_h2 <- gsim:::.gsim_packed_unpack(reloaded$h2)
  testthat::expect_identical(loaded_h1, raw_pedigree$h1)
  testthat::expect_identical(loaded_h2, raw_pedigree$h2)
  testthat::expect_equal(.hap_mendelian_errors(
    loaded_h1, loaded_h2, raw_pedigree$pedigree_alignment), 0L)

  bed <- gsim:::.gsim_bed_sink_create(
    backend, file.path(root, "reloaded.bed"), rownames(loaded_h1))
  gsim:::.gsim_bed_sink_append(bed, "ped", reloaded$h1, reloaded$h2, pvars)
  bed_manifest <- gsim:::.gsim_bed_sink_finalize(bed)
  decoded <- gsim:::.gsim_packed_bed_read_all(
    backend, bed_manifest$path, nrow(loaded_h1), ncol(loaded_h1),
    rownames(loaded_h1), colnames(loaded_h1))
  expected_genotypes <- matrix(as.integer(raw_pedigree$genotypes),
                               nrow(raw_pedigree$genotypes),
                               dimnames = dimnames(raw_pedigree$genotypes))
  testthat::expect_identical(decoded, expected_genotypes)
  testthat::expect_identical(decoded,
    matrix(as.integer(loaded_h1) + as.integer(loaded_h2), nrow(loaded_h1),
           dimnames = dimnames(loaded_h1)))
  gsim:::.gsim_hap_dataset_close(pedigree_reader)
})
