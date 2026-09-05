.href_backend <- function() {
  gsim:::.gsim_packed_backend()
}

.href_metadata_backend <- function() {
  gsim:::.gsim_metadata_backend()
}

.href_variants <- function(label, ids, position) {
  data.frame(
    chromosome = rep.int(label, length(ids)), variant_id = ids,
    genetic_position_cm = position,
    base_pair_position = seq_len(length(ids)),
    alt = rep.int("C", length(ids)), ref = rep.int("A", length(ids)),
    bit1_allele = rep.int("C", length(ids)),
    bit0_allele = rep.int("A", length(ids)), stringsAsFactors = FALSE)
}

.href_panel <- function(markers = c(`Z-2` = 13L, `01` = 7L, chrA = 9L)) {
  donors <- paste0("d", seq_len(4L))
  panels <- lapply(seq_along(markers), function(block) {
    m <- markers[[block]]
    ids <- paste0(names(markers)[[block]], "_v", seq_len(m))
    h1 <- outer(seq_along(donors), seq_len(m), function(i, j) {
      (i * 3L + j * 5L + block) %% 2L
    })
    h2 <- outer(seq_along(donors), seq_len(m), function(i, j) {
      (i * 7L + j * 3L + j %/% 2L + block + 1L) %% 2L
    })
    dimnames(h1) <- dimnames(h2) <- list(donors, ids)
    position <- seq(0, 2 + block / 10, length.out = m)
    list(h1 = h1, h2 = h2, ids = ids, position = position,
         mutation = rep(c(0.25, 2, 20, 1e9), length.out = m))
  })
  names(panels) <- names(markers)
  panels
}

.href_write <- function(backend, metadata_backend, prefix, panels) {
  donors <- rownames(panels[[1L]]$h1)
  dataset <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, prefix,
    gsim:::.gsim_plink_sample_metadata(donors))
  for (label in names(panels)) {
    panel <- panels[[label]]
    h1 <- gsim:::.gsim_packed_pack(backend, panel$h1)
    h2 <- gsim:::.gsim_packed_pack(backend, panel$h2)
    gsim:::.gsim_hap_dataset_append(
      dataset, label, h1, h2,
      .href_variants(label, panel$ids, panel$position))
    gsim:::.gsim_packed_close(h1)
    gsim:::.gsim_packed_close(h2)
  }
  gsim:::.gsim_hap_dataset_finalize(dataset)
}

.href_args <- function(panel, label, n = 8L, seed = 717,
                       individual_offset = 0L) {
  list(
    donor_population = c(d1 = "A", d2 = "A", d3 = "B", d4 = "B"),
    ancestry_weights = c(A = 0.35, B = 0.65),
    N = c(A = 2, B = 2), Ne = c(A = 4, B = 7),
    rho = c(A = 0.8, B = 1.3),
    genetic_position = stats::setNames(panel$position, panel$ids),
    mutation_age = stats::setNames(panel$mutation, panel$ids),
    n = n, seed = seed, chromosome = label,
    individual_offset = individual_offset)
}

.href_raw <- function(panel, label, args, return_genotypes = TRUE,
                      return_segments = TRUE) {
  raw_args <- args
  raw_args$donor_population <- unname(raw_args$donor_population)
  raw_args$genetic_position <- unname(raw_args$genetic_position)
  raw_args$mutation_age <- unname(raw_args$mutation_age)
  raw_args$chromosome <- rep.int(label, ncol(panel$h1))
  do.call(gsim:::.gsim_hapnest_founders, c(
    list(reference_haplotypes_h1 = panel$h1,
         reference_haplotypes_h2 = panel$h2), raw_args,
    list(return_genotypes = return_genotypes,
         return_segments = return_segments)))
}

.href_simulate <- function(reader, label, panel, args,
                           return_genotypes = TRUE, return_segments = TRUE) {
  do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome, c(
    list(dataset = reader), args,
    list(return_genotypes = return_genotypes,
         return_segments = return_segments)))
}

.href_mendelian_errors <- function(h1, h2, alignment) {
  errors <- 0L
  for (i in which(!alignment$founder)) {
    sire <- match(alignment$sire[[i]], alignment$animal)
    dam <- match(alignment$dam[[i]], alignment$animal)
    errors <- errors + sum(!(h1[i, ] == h1[sire, ] | h1[i, ] == h2[sire, ]))
    errors <- errors + sum(!(h2[i, ] == h1[dam, ] | h2[i, ] == h2[dam, ]))
  }
  errors
}

testthat::test_that("HAP-loaded packed donors exactly match both founder oracles", {
  backend <- .href_backend()
  metadata_backend <- .href_metadata_backend()
  root <- tempfile("hap-reference-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  panels <- .href_panel()
  prefix <- file.path(root, "reference panel")
  .href_write(backend, metadata_backend, prefix, panels)
  reader <- gsim:::.gsim_hap_dataset_open(backend, metadata_backend, prefix)
  on.exit(gsim:::.gsim_hap_dataset_close(reader), add = TRUE)

  panel <- panels[["Z-2"]]
  args <- .href_args(panel, "Z-2")
  raw <- .href_raw(panel, "Z-2", args)
  legacy <- do.call(gsim:::.gsim_hapnest_founders_packed_chromosome, c(
    list(backend = backend, reference_haplotypes_h1 = panel$h1,
         reference_haplotypes_h2 = panel$h2), args,
    list(return_genotypes = TRUE, return_segments = TRUE)))
  loaded <- .href_simulate(reader, "Z-2", panel, args)
  permuted_args <- args
  permuted_args$donor_population <- rev(permuted_args$donor_population)
  permuted_args$genetic_position <- rev(permuted_args$genetic_position)
  permuted_args$mutation_age <- rev(permuted_args$mutation_age)
  permuted <- .href_simulate(reader, "Z-2", panel, permuted_args)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(loaded$h1), raw$h1)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(loaded$h2), raw$h2)
  testthat::expect_identical(loaded$genotypes, raw$genotypes)
  testthat::expect_identical(loaded$segments, raw$segments)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(permuted$h1), raw$h1)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(permuted$h2), raw$h2)
  testthat::expect_identical(permuted$segments, raw$segments)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(loaded$h1),
                             gsim:::.gsim_packed_unpack(legacy$h1))
  testthat::expect_identical(gsim:::.gsim_packed_unpack(loaded$h2),
                             gsim:::.gsim_packed_unpack(legacy$h2))
  testthat::expect_identical(loaded$segments, legacy$segments)
  testthat::expect_identical(
    as.integer(loaded$genotypes),
    as.integer(gsim:::.gsim_packed_unpack(loaded$h1)) +
      as.integer(gsim:::.gsim_packed_unpack(loaded$h2)))
  testthat::expect_identical(loaded$reference_alignment$sample_ids,
                             rownames(panel$h1))
  testthat::expect_identical(loaded$reference_alignment$variant_ids, panel$ids)

  generated_prefix <- file.path(root, "generated founders")
  generated_dataset <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, generated_prefix,
    gsim:::.gsim_plink_sample_metadata(loaded$sample_ids))
  gsim:::.gsim_hap_dataset_append(
    generated_dataset, "Z-2", loaded$h1, loaded$h2,
    .href_variants("Z-2", panel$ids, panel$position))
  gsim:::.gsim_hap_dataset_finalize(generated_dataset)
  generated_reader <- gsim:::.gsim_hap_dataset_open(
    backend, metadata_backend, generated_prefix)
  reloaded <- gsim:::.gsim_hap_dataset_load_chromosome(
    generated_reader, "Z-2")
  testthat::expect_identical(gsim:::.gsim_packed_unpack(reloaded$h1), raw$h1)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(reloaded$h2), raw$h2)
  gsim:::.gsim_hap_dataset_close(generated_reader)
  bed <- gsim:::.gsim_bed_sink_create(
    backend, file.path(root, "generated founders.bed"), loaded$sample_ids)
  gsim:::.gsim_bed_sink_append(
    bed, "Z-2", loaded$h1, loaded$h2, loaded$variant_ids)
  bed_manifest <- gsim:::.gsim_bed_sink_finalize(bed)
  decoded <- gsim:::.gsim_packed_bed_read_all(
    backend, bed_manifest$path, length(loaded$sample_ids),
    length(loaded$variant_ids), loaded$sample_ids, loaded$variant_ids)
  testthat::expect_identical(
    decoded, matrix(as.integer(raw$genotypes), nrow(raw$genotypes),
                    dimnames = dimnames(raw$genotypes)))
})

testthat::test_that("HAP-loaded copying retains strict phase specificity", {
  backend <- .href_backend()
  metadata_backend <- .href_metadata_backend()
  root <- tempfile("hap-phase-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  h1 <- matrix(0, 2L, 11L,
               dimnames = list(c("p1", "p2"), paste0("q", 1:11)))
  h2 <- matrix(1, 2L, 11L, dimnames = dimnames(h1))
  panel <- list(h1 = h1, h2 = h2, ids = colnames(h1),
                position = seq(0, 1, length.out = 11L),
                mutation = rep(1e12, 11L))
  prefix <- file.path(root, "phase")
  .href_write(backend, metadata_backend, prefix, list(asym = panel))
  reader <- gsim:::.gsim_hap_dataset_open(backend, metadata_backend, prefix)
  on.exit(gsim:::.gsim_hap_dataset_close(reader), add = TRUE)
  args <- list(
    donor_population = c(p2 = "P", p1 = "P"), ancestry_weights = c(P = 1),
    N = c(P = 2), Ne = c(P = 2), rho = c(P = 1),
    genetic_position = stats::setNames(panel$position, panel$ids),
    mutation_age = stats::setNames(panel$mutation, panel$ids),
    n = 6L, seed = 19L, chromosome = "asym")
  out <- .href_simulate(reader, "asym", panel, args)
  testthat::expect_true(all(gsim:::.gsim_packed_unpack(out$h1) == 0))
  testthat::expect_true(all(gsim:::.gsim_packed_unpack(out$h2) == 1))
  testthat::expect_identical(out$settings$donor_phase, "hapnest")
})

testthat::test_that("chromosome, batching, options, and handle lifetime are invariant", {
  backend <- .href_backend()
  metadata_backend <- .href_metadata_backend()
  root <- tempfile("hap-invariance-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  panels <- .href_panel()
  multi <- file.path(root, "multi")
  .href_write(backend, metadata_backend, multi, panels)
  reader <- gsim:::.gsim_hap_dataset_open(backend, metadata_backend, multi)
  forward <- lapply(names(panels), function(label) {
    .href_simulate(reader, label, panels[[label]], .href_args(panels[[label]], label))
  })
  names(forward) <- names(panels)
  reverse <- lapply(rev(names(panels)), function(label) {
    .href_simulate(reader, label, panels[[label]], .href_args(panels[[label]], label))
  })
  names(reverse) <- rev(names(panels))
  for (label in names(panels)) {
    testthat::expect_identical(gsim:::.gsim_packed_unpack(forward[[label]]$h1),
                               gsim:::.gsim_packed_unpack(reverse[[label]]$h1))
    testthat::expect_identical(gsim:::.gsim_packed_unpack(forward[[label]]$h2),
                               gsim:::.gsim_packed_unpack(reverse[[label]]$h2))
    testthat::expect_identical(forward[[label]]$segments,
                               reverse[[label]]$segments)
  }
  standalone <- file.path(root, "standalone")
  .href_write(backend, metadata_backend, standalone, panels["01"])
  one_reader <- gsim:::.gsim_hap_dataset_open(backend, metadata_backend, standalone)
  alone <- .href_simulate(one_reader, "01", panels[["01"]],
                          .href_args(panels[["01"]], "01"))
  gsim:::.gsim_hap_dataset_close(one_reader)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(alone$h1),
                             gsim:::.gsim_packed_unpack(forward[["01"]]$h1))
  testthat::expect_identical(alone$segments, forward[["01"]]$segments)

  panel <- panels[["Z-2"]]
  full <- forward[["Z-2"]]
  first <- .href_simulate(reader, "Z-2", panel,
                          .href_args(panel, "Z-2", n = 3L))
  second <- .href_simulate(reader, "Z-2", panel,
                           .href_args(panel, "Z-2", n = 5L,
                                      individual_offset = 3L))
  testthat::expect_identical(
    rbind(gsim:::.gsim_packed_unpack(first$h1),
          gsim:::.gsim_packed_unpack(second$h1)),
    gsim:::.gsim_packed_unpack(full$h1))
  segments <- rbind(first$segments, second$segments)
  rownames(segments) <- NULL
  expected_segments <- full$segments
  rownames(expected_segments) <- NULL
  testthat::expect_identical(segments, expected_segments)
  without_outputs <- .href_simulate(
    reader, "Z-2", panel, .href_args(panel, "Z-2"),
    return_genotypes = FALSE, return_segments = FALSE)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(without_outputs$h1),
                             gsim:::.gsim_packed_unpack(full$h1))
  testthat::expect_null(without_outputs$genotypes)
  testthat::expect_null(without_outputs$segments)
  changed_seed <- .href_simulate(
    reader, "Z-2", panel, .href_args(panel, "Z-2", seed = 718L),
    return_genotypes = FALSE)
  testthat::expect_false(identical(
    gsim:::.gsim_packed_unpack(changed_seed$h1),
    gsim:::.gsim_packed_unpack(full$h1)))
  set.seed(41)
  expected_rng <- runif(5)
  set.seed(41)
  invisible(.href_simulate(reader, "01", panels[["01"]],
                            .href_args(panels[["01"]], "01")))
  testthat::expect_identical(runif(5), expected_rng)

  owned <- gsim:::.gsim_hap_dataset_load_chromosome(reader, "chrA")
  gsim:::.gsim_hap_dataset_close(reader)
  after_close <- do.call(gsim:::.gsim_hapnest_founders_packed_reference_chromosome,
    c(list(backend = backend, reference_h1 = owned$h1, reference_h2 = owned$h2),
      .href_args(panels$chrA, "chrA")))
  testthat::expect_identical(gsim:::.gsim_packed_unpack(after_close$h1),
                             gsim:::.gsim_packed_unpack(forward$chrA$h1))
  gsim:::.gsim_packed_close(owned$h1)
  gsim:::.gsim_packed_close(owned$h2)
  testthat::expect_error(gsim:::.gsim_packed_info(owned$h1), "released")
  testthat::expect_error(gsim:::.gsim_packed_close(owned$h1), "released")
})

testthat::test_that("reference alignment validation rejects ambiguity", {
  backend <- .href_backend()
  metadata_backend <- .href_metadata_backend()
  root <- tempfile("hap-validation-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  panels <- .href_panel(c(chrA = 9L))
  panel <- panels$chrA
  prefix <- file.path(root, "reference")
  .href_write(backend, metadata_backend, prefix, panels)
  reader <- gsim:::.gsim_hap_dataset_open(backend, metadata_backend, prefix)
  on.exit(gsim:::.gsim_hap_dataset_close(reader), add = TRUE)
  args <- .href_args(panel, "chrA")
  bad <- args
  bad$donor_population <- bad$donor_population[-1L]
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
                                 c(list(dataset = reader), bad)), "missing: d1")
  bad <- args
  names(bad$donor_population)[2L] <- "d1"
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
                                 c(list(dataset = reader), bad)), "exactly match")
  bad <- args
  bad$donor_population <- c(bad$donor_population, extra = "A")
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
                                 c(list(dataset = reader), bad)), "extra: extra")
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
    c(list(dataset = reader), within(args, chromosome <- "absent"))), "identify")
  bad <- args
  bad$genetic_position[[2L]] <- bad$genetic_position[[2L]] + 0.01
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
                                 c(list(dataset = reader), bad)), "serialized BIM")
  bad <- args
  bad$mutation_age <- bad$mutation_age[-1L]
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
                                 c(list(dataset = reader), bad)), "missing")
  bad <- args
  bad$ancestry_weights <- c(A = 0, B = 0)
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
                                 c(list(dataset = reader), bad)), "positive")
  bad <- args
  bad$N <- c(A = 1, B = 2)
  testthat::expect_error(do.call(gsim:::.gsim_hapnest_founders_from_hap_chromosome,
                                 c(list(dataset = reader), bad)), "N must equal")

  loaded <- gsim:::.gsim_hap_dataset_load_chromosome(reader, "chrA")
  wrong <- matrix(0, 4L, 8L,
                  dimnames = list(rownames(panel$h1), panel$ids[-1L]))
  wrong_h2 <- gsim:::.gsim_packed_pack(backend, wrong)
  testthat::expect_error(do.call(
    gsim:::.gsim_hapnest_founders_packed_reference_chromosome,
    c(list(backend = backend, reference_h1 = loaded$h1, reference_h2 = wrong_h2),
      args)), "dimensions")
  gsim:::.gsim_packed_close(wrong_h2)
  gsim:::.gsim_packed_close(loaded$h1)
  gsim:::.gsim_packed_close(loaded$h2)
})

testthat::test_that("generated packed founders flow directly through pedigree, HAP, and BED", {
  backend <- .href_backend()
  metadata_backend <- .href_metadata_backend()
  root <- tempfile("hap-downstream-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  panels <- .href_panel(c(chrP = 25L))
  panel <- panels$chrP
  prefix <- file.path(root, "reference")
  .href_write(backend, metadata_backend, prefix, panels)
  reader <- gsim:::.gsim_hap_dataset_open(backend, metadata_backend, prefix)
  args <- .href_args(panel, "chrP", n = 5L, seed = 919)
  founders <- .href_simulate(reader, "chrP", panel, args,
                             return_genotypes = FALSE)
  raw_founders <- .href_raw(panel, "chrP", args, return_genotypes = FALSE)
  gsim:::.gsim_hap_dataset_close(reader)

  tab <- data.frame(
    animal = c("syn1", "syn2", "syn3", "syn4", "C", "FS", "PHS",
               "MHS", "G", "syn5"),
    sire = c(rep(NA, 4), "syn1", "syn1", "syn1", "syn4", "C", NA),
    dam = c(rep(NA, 4), "syn2", "syn2", "syn3", "syn2", "syn3", NA),
    stringsAsFactors = FALSE)
  external <- c(10, 7, 3, 1, 9, 2, 5, 8, 4, 6)
  pedigree <- structure(list(
    pedigree = tab[external, , drop = FALSE], canonical_order = tab$animal,
    external_order = tab$animal[external]), class = "gsim_pedigree")
  meiosis_args <- list(
    pedigree = pedigree, chromosome = rep.int("chrP", ncol(panel$h1)),
    genetic_position = seq(0, 8, length.out = ncol(panel$h1)),
    seed = 303L, return_genotypes = TRUE, return_crossovers = TRUE)
  packed <- do.call(gsim:::.gsim_pedigree_genotypes_packed_chromosome, c(
    list(backend = backend,
         founder_haplotypes = list(h1 = founders$h1, h2 = founders$h2)),
    meiosis_args))
  raw <- do.call(gsim:::.gsim_pedigree_genotypes, c(
    list(founder_haplotypes = list(h1 = raw_founders$h1,
                                   h2 = raw_founders$h2)), meiosis_args))
  packed_h1 <- gsim:::.gsim_packed_unpack(packed$h1)
  packed_h2 <- gsim:::.gsim_packed_unpack(packed$h2)
  testthat::expect_identical(packed_h1, raw$h1)
  testthat::expect_identical(packed_h2, raw$h2)
  testthat::expect_identical(packed$genotypes, raw$genotypes)
  testthat::expect_identical(packed$crossover_audit, raw$crossover_audit)
  testthat::expect_equal(.href_mendelian_errors(
    packed_h1, packed_h2, packed$pedigree_alignment), 0L)

  generated_prefix <- file.path(root, "generated pedigree")
  dataset <- gsim:::.gsim_hap_dataset_create(
    backend, metadata_backend, generated_prefix,
    gsim:::.gsim_plink_pedigree_metadata(pedigree))
  gsim:::.gsim_hap_dataset_append(
    dataset, "chrP", packed$h1, packed$h2,
    .href_variants("chrP", panel$ids, seq(0, 8, length.out = length(panel$ids))))
  gsim:::.gsim_hap_dataset_finalize(dataset)
  reread <- gsim:::.gsim_hap_dataset_open(
    backend, metadata_backend, generated_prefix)
  loaded <- gsim:::.gsim_hap_dataset_load_chromosome(reread, "chrP")
  testthat::expect_identical(gsim:::.gsim_packed_unpack(loaded$h1), raw$h1)
  testthat::expect_identical(gsim:::.gsim_packed_unpack(loaded$h2), raw$h2)
  gsim:::.gsim_hap_dataset_close(reread)

  bed <- gsim:::.gsim_bed_sink_create(
    backend, file.path(root, "generated.bed"), packed$sample_ids)
  gsim:::.gsim_bed_sink_append(
    bed, "chrP", packed$h1, packed$h2, packed$variant_ids)
  bed_manifest <- gsim:::.gsim_bed_sink_finalize(bed)
  decoded <- gsim:::.gsim_packed_bed_read_all(
    backend, bed_manifest$path, length(packed$sample_ids),
    length(packed$variant_ids), packed$sample_ids, packed$variant_ids)
  testthat::expect_identical(
    decoded,
    matrix(as.integer(raw$genotypes), nrow(raw$genotypes),
           dimnames = dimnames(raw$genotypes)))
})
