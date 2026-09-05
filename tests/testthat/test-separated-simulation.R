.separated_write_vcf <- function(path, samples, chromosomes = c("B", "A"),
                                 markers = 16L) {
  header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", samples), collapse = "\t")
  patterns <- c("0|0", "0|1", "1|0", "1|1")
  records <- unlist(lapply(chromosomes, function(chromosome) {
    vapply(seq_len(markers), function(marker) {
      calls <- patterns[((seq_along(samples) + marker - 2L) %% 4L) + 1L]
      paste(c(chromosome, marker * 10L,
              paste0(tolower(chromosome), marker), "A", "G", ".", ".",
              ".", "GT", calls), collapse = "\t")
    }, character(1L))
  }), use.names = FALSE)
  writeBin(charToRaw(paste(c("##fileformat=VCFv4.2", header, records, ""),
                           collapse = "\n")), path)
}

.separated_pedigree <- function(founder_ids) {
  children <- c("FS1", "FS2", "PHS", "MHS", "G")
  tab <- data.frame(
    animal = c(founder_ids, children),
    sire = c(rep(NA_character_, length(founder_ids)),
             founder_ids[[1L]], founder_ids[[1L]], founder_ids[[1L]],
             founder_ids[[4L]], "FS1"),
    dam = c(rep(NA_character_, length(founder_ids)),
            founder_ids[[2L]], founder_ids[[2L]], founder_ids[[3L]],
            founder_ids[[2L]], founder_ids[[5L]]),
    sex = "U",
    generation = c(rep(1L, length(founder_ids)), 2L, 2L, 2L, 2L, 3L),
    cohort = c(rep(1L, length(founder_ids)), 2L, 2L, 2L, 2L, 3L),
    phenotyped = TRUE, stringsAsFactors = FALSE
  )
  canonical <- tab$animal
  external <- rev(canonical)
  structure(list(
    pedigree = tab[match(external, canonical), , drop = FALSE],
    canonical_order = canonical, external_order = external,
    mapping = data.frame(
      animal = canonical, canonical_index = seq_along(canonical),
      external_index = match(canonical, external)
    ), settings = list(), diagnostics = list(), checksums = list()
  ), class = "gsim_pedigree")
}

.separated_read_phases <- function(prefix) {
  packed <- NULL
  metadata <- NULL
  reader <- gsim:::.gsim_hap_dataset_open(prefix)
  on.exit(gsim:::.gsim_hap_dataset_close(reader), add = TRUE)
  setNames(lapply(reader$chromosome, function(chromosome) {
    phase <- gsim:::.gsim_hap_dataset_load_chromosome(reader, chromosome)
    on.exit({
      gsim:::.gsim_packed_close(phase$h1)
      gsim:::.gsim_packed_close(phase$h2)
    }, add = TRUE)
    list(h1 = gsim:::.gsim_packed_unpack(phase$h1),
         h2 = gsim:::.gsim_packed_unpack(phase$h2))
  }), reader$chromosome)
}

testthat::test_that("separated founder batches and threads are exactly invariant", {
  root <- tempfile("separated-simulation-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  samples <- paste0("d", 1:4)
  vcf <- file.path(root, "reference.vcf")
  .separated_write_vcf(vcf, samples)
  variant_ids <- unlist(lapply(c("B", "A"), function(x) {
    paste0(tolower(x), seq_len(16L))
  }), use.names = FALSE)
  map <- data.frame(
    chromosome = rep(c("B", "A"), each = 16L),
    variant_id = variant_ids,
    genetic_position_cm = rep(seq(0, 150, length.out = 16L), 2L),
    stringsAsFactors = FALSE
  )
  reference <- gsim_import_vcf(vcf, map, file.path(root, "reference"))
  populations <- stats::setNames(c("P1", "P1", "P2", "P2"), samples)
  mutation_age <- stats::setNames(rep(1e9, length(variant_ids)), variant_ids)
  founder_arguments <- list(
    reference = reference, populations = populations,
    ancestry_weights = c(P1 = 0.4, P2 = 0.6),
    mutation_age = mutation_age, N = c(P1 = 2, P2 = 2),
    Ne = c(P1 = 4, P2 = 6), rho = c(P1 = 0.02, P2 = 0.03),
    seed = 717, output = file.path(root, "invalid-base")
  )
  testthat::expect_error(
    do.call(gsim_simulate_founders, founder_arguments), "Exactly one"
  )
  testthat::expect_error(
    do.call(gsim_simulate_founders,
            c(founder_arguments, list(n = 2L, founder_ids = c("x", "y")))),
    "Exactly one"
  )
  testthat::expect_error(
    do.call(gsim_simulate_founders,
            c(founder_arguments, list(founder_ids = c("x", "x")))),
    "unique"
  )
  testthat::expect_error(
    do.call(gsim_simulate_founders,
            c(founder_arguments, list(founder_ids = 1:2))),
    "character"
  )
  founder_ids <- sprintf("F%03d", seq_len(513L))
  configurations <- data.frame(
    batch_size = c(10000L, 70L, 64L, 129L),
    threads = c(1L, 2L, 4L, 8L)
  )
  before <- { set.seed(741); .Random.seed }
  bases <- lapply(seq_len(nrow(configurations)), function(i) {
    gsim:::.gsim_simulate_founders_impl(
      reference, founder_ids = founder_ids, populations = populations,
      ancestry_weights = c(P1 = 0.4, P2 = 0.6), mutation_age = mutation_age,
      N = c(P1 = 2, P2 = 2), Ne = c(P1 = 4, P2 = 6),
      rho = c(P1 = 0.02, P2 = 0.03), seed = 717,
      output = file.path(root, paste0("base-", i)),
      batch_size = configurations$batch_size[[i]],
      threads = configurations$threads[[i]], return_segments = TRUE
    )
  })
  testthat::expect_identical(.Random.seed, before)
  hap_bytes <- lapply(bases, function(x) readBin(x$paths[["hap"]], "raw",
                                                 file.info(x$paths[["hap"]])$size))
  for (i in 2:length(bases)) {
    testthat::expect_identical(hap_bytes[[i]], hap_bytes[[1L]])
    testthat::expect_identical(bases[[i]]$sample_ids, founder_ids)
    testthat::expect_identical(bases[[i]]$segment_audit,
                               bases[[1L]]$segment_audit)
    testthat::expect_identical(.separated_read_phases(bases[[i]]$prefix),
                               .separated_read_phases(bases[[1L]]$prefix))
  }
  testthat::expect_equal(bases[[2L]]$simulation$actual_batch_size, 128L)
  testthat::expect_equal(bases[[3L]]$simulation$actual_batch_size, 64L)
  testthat::expect_lt(
    bases[[3L]]$simulation$peak_chromosome_biological_payload_bytes,
    2 * 8 * 16 * ceiling(length(founder_ids) / 64) +
      2 * 8 * 16 * ceiling(length(samples) / 64)
  )

  pedigree <- .separated_pedigree(founder_ids)
  base_hash_before <- unname(tools::md5sum(bases[[1L]]$paths[["hap"]]))
  pedigree_hap <- gsim:::.gsim_simulate_pedigree_impl(
    bases[[1L]], pedigree, 991, file.path(root, "pedigree-hap"), "hap",
    return_crossovers = TRUE
  )
  pedigree_bed <- gsim_simulate_pedigree(
    bases[[4L]], pedigree, 991, file.path(root, "pedigree-bed"), "bed"
  )
  second_pedigree <- gsim_simulate_pedigree(
    bases[[1L]], pedigree, 992, file.path(root, "pedigree-second"), "hap"
  )
  testthat::expect_identical(
    unname(tools::md5sum(bases[[1L]]$paths[["hap"]])), base_hash_before
  )
  first_phases <- .separated_read_phases(file.path(root, "pedigree-hap"))
  second_phases <- .separated_read_phases(file.path(root, "pedigree-second"))
  testthat::expect_false(identical(first_phases, second_phases))
  packed <- NULL
  decoded <- gsim:::.gsim_packed_bed_read_all(
    pedigree_bed$paths[["bed"]], length(pedigree$canonical_order),
    length(variant_ids)
  )
  expected <- do.call(cbind, lapply(first_phases, function(x) {
    matrix(as.integer(x$h1) + as.integer(x$h2), nrow = nrow(x$h1))
  }))
  testthat::expect_identical(unname(decoded), unname(expected))
  inconsistency <- 0L
  alignment <- pedigree$pedigree[
    match(pedigree$canonical_order, pedigree$pedigree$animal), , drop = FALSE
  ]
  for (phase in first_phases) {
    for (i in which(!is.na(alignment$sire))) {
      sire <- match(alignment$sire[[i]], pedigree$canonical_order)
      dam <- match(alignment$dam[[i]], pedigree$canonical_order)
      inconsistency <- inconsistency +
        sum(!(phase$h1[i, ] == phase$h1[sire, ] |
              phase$h1[i, ] == phase$h2[sire, ])) +
        sum(!(phase$h2[i, ] == phase$h1[dam, ] |
              phase$h2[i, ] == phase$h2[dam, ]))
    }
  }
  testthat::expect_equal(inconsistency, 0L)
  testthat::expect_error(
    gsim_simulate_pedigree(
      bases[[1L]], .separated_pedigree(c(founder_ids[-1L], "missing")),
      991, file.path(root, "bad-map")
    ), "must exactly match"
  )
  testthat::expect_identical(pedigree_hap$sample_ids,
                             pedigree$canonical_order)
})
