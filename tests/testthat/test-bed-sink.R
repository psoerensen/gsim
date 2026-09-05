.bed_backend <- function() {
  gsim:::.gsim_packed_backend()
}

.bed_expected <- function(genotypes) {
  stopifnot(is.matrix(genotypes))
  genotype_values <- matrix(as.integer(genotypes), nrow(genotypes), ncol(genotypes))
  stopifnot(all(genotype_values %in% 0:2))
  individuals <- nrow(genotypes)
  markers <- ncol(genotypes)
  bytes_per_marker <- ceiling(individuals / 4)
  output <- as.raw(c(0x6c, 0x1b, 0x01, rep.int(0L, markers * bytes_per_marker)))
  for (marker in seq_len(markers)) {
    for (individual in seq_len(individuals)) {
      code <- c(3L, 2L, 0L)[genotype_values[individual, marker] + 1L]
      at <- 3L + (marker - 1L) * bytes_per_marker +
        ((individual - 1L) %/% 4L) + 1L
      shift <- 2L * ((individual - 1L) %% 4L)
      output[[at]] <- as.raw(bitwOr(as.integer(output[[at]]),
                                   bitwShiftL(code, shift)))
    }
  }
  output
}

.bed_founder_fixture <- function(backend, chromosome = "chrZ",
                                 markers = 11L, n = 7L, seed = 1905) {
  marker <- seq_len(markers)
  h1 <- outer(seq_len(4L), marker,
              function(i, j) (i * 3L + j + j %/% 2L) %% 2L)
  h2 <- outer(seq_len(4L), marker,
              function(i, j) (i + j * 5L + j %/% 3L + 1L) %% 2L)
  dimnames(h1) <- dimnames(h2) <- list(
    paste0("donor", seq_len(4L)), paste0(chromosome, "_v", marker)
  )
  args <- list(
    reference_haplotypes_h1 = h1,
    reference_haplotypes_h2 = h2,
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

.bed_pedigree_fixture <- function(backend, chromosome = "chrP",
                                  markers = 19L, seed = 905) {
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
  founder_ids <- c("S", "D", "D2", "S2", "L")
  variants <- paste0(chromosome, "_v", seq_len(markers))
  h1 <- outer(seq_along(founder_ids), seq_len(markers),
              function(i, j) (i + j + j %/% 4L) %% 2L)
  h2 <- outer(seq_along(founder_ids), seq_len(markers),
              function(i, j) (i * 3L + j * 5L + j %/% 3L + 1L) %% 2L)
  dimnames(h1) <- dimnames(h2) <- list(founder_ids, variants)
  args <- list(
    pedigree = pedigree, founder_haplotypes = list(h1 = h1, h2 = h2),
    chromosome = rep(chromosome, markers),
    genetic_position = seq(0, 7, length.out = markers),
    seed = seed, return_crossovers = TRUE
  )
  list(
    raw = do.call(gsim:::.gsim_pedigree_genotypes,
                  c(args, list(return_genotypes = TRUE))),
    packed = do.call(gsim:::.gsim_pedigree_genotypes_packed_chromosome,
                     c(list(backend = backend), args,
                       list(return_genotypes = FALSE)))
  )
}

.write_one_bed <- function(backend, path, packed, chromosome,
                           buffer_variants = 2L, overwrite = FALSE) {
  sink <- gsim:::.gsim_bed_sink_create(
    backend, path, packed$sample_ids, overwrite = overwrite,
    buffer_variants = buffer_variants,
    provenance = list(seed = packed$settings$seed, model = packed$settings$model)
  )
  gsim:::.gsim_bed_sink_append(
    sink, chromosome, packed$h1, packed$h2, packed$variant_ids
  )
  gsim:::.gsim_bed_sink_finalize(sink)
}

testthat::test_that("R sink preserves hand-calculated SNP-major bytes", {
  backend <- .bed_backend()
  samples <- paste0("i", 1:4)
  variants <- "v1"
  h1 <- matrix(as.raw(c(0, 0, 1, 0)), 4L, 1L,
               dimnames = list(samples, variants))
  h2 <- matrix(as.raw(c(0, 1, 1, 0)), 4L, 1L,
               dimnames = list(samples, variants))
  path <- tempfile(fileext = ".bed")
  sink <- gsim:::.gsim_bed_sink_create(backend, path, samples,
                                       buffer_variants = 1L)
  gsim:::.gsim_bed_sink_append(
    sink, "chr1", gsim:::.gsim_packed_pack(backend, h1),
    gsim:::.gsim_packed_pack(backend, h2), variants
  )
  manifest <- gsim:::.gsim_bed_sink_finalize(sink)
  testthat::expect_identical(readBin(path, "raw", n = 99L),
                             as.raw(c(0x6c, 0x1b, 0x01, 0xcb)))
  testthat::expect_identical(
    gsim:::.gsim_packed_bed_read_all(backend, path, 4L, 1L,
                                    samples, variants),
    matrix(c(0L, 1L, 2L, 0L), 4L, 1L,
           dimnames = list(samples, variants))
  )
  testthat::expect_identical(manifest$bytes_written, 4)
  testthat::expect_identical(manifest$expected_bytes, 4)
  testthat::expect_identical(manifest$implementation$engine,
                             "gsim private native backend")
  testthat::expect_match(manifest$implementation$origin, "089bf1e")
  testthat::expect_identical(manifest$chromosome_order, "chr1")
  set.seed(911)
  expected_rng <- runif(4)
  set.seed(911)
  invisible(gsim:::.gsim_packed_bed_read_all(backend, path, 4L, 1L))
  testthat::expect_identical(runif(4), expected_rng)
  unlink(path)
})

testthat::test_that("founder BED decode exactly matches both qualified oracles", {
  backend <- .bed_backend()
  fixture <- .bed_founder_fixture(backend)
  path <- tempfile(fileext = ".bed")
  manifest <- .write_one_bed(backend, path, fixture$packed, "chrZ", 1L)
  decoded <- gsim:::.gsim_packed_bed_read_all(
    backend, path, length(fixture$packed$sample_ids),
    length(fixture$packed$variant_ids), fixture$packed$sample_ids,
    fixture$packed$variant_ids
  )
  bounded <- gsim:::.gsim_packed_decode_genotypes(
    fixture$packed$h1, fixture$packed$h2
  )
  testthat::expect_identical(decoded, matrix(as.integer(fixture$raw$genotypes),
                                             nrow(fixture$raw$genotypes),
                                             dimnames = dimnames(fixture$raw$genotypes)))
  testthat::expect_identical(decoded,
                             matrix(as.integer(bounded), nrow(bounded),
                                    dimnames = dimnames(bounded)))
  unpacked_sum <- matrix(
    as.integer(gsim:::.gsim_packed_unpack(fixture$packed$h1)) +
      as.integer(gsim:::.gsim_packed_unpack(fixture$packed$h2)),
    nrow(decoded), dimnames = dimnames(decoded)
  )
  testthat::expect_identical(decoded, unpacked_sum)
  testthat::expect_identical(sum(decoded == -9L), 0L)
  testthat::expect_identical(readBin(path, "raw", n = file.info(path)$size),
                             .bed_expected(fixture$raw$genotypes))
  testthat::expect_identical(manifest$variant_ids,
                             colnames(fixture$raw$genotypes))
  testthat::expect_identical(manifest$sample_ids,
                             rownames(fixture$raw$genotypes))
  unlink(path)

  h1 <- matrix(as.raw(0), 7L, 3L,
               dimnames = list(paste0("p", 1:7), paste0("a", 1:3)))
  h2 <- matrix(as.raw(1), 7L, 3L, dimnames = dimnames(h1))
  asymmetric <- list(
    h1 = gsim:::.gsim_packed_pack(backend, h1),
    h2 = gsim:::.gsim_packed_pack(backend, h2),
    sample_ids = rownames(h1), variant_ids = colnames(h1),
    settings = list(seed = 1, model = "phase-asymmetric fixture")
  )
  path <- tempfile(fileext = ".bed")
  .write_one_bed(backend, path, asymmetric, "asym", 2L)
  testthat::expect_identical(readBin(path, "raw", n = 99L),
                             as.raw(c(0x6c, 0x1b, 0x01,
                                      rep(c(0xaa, 0x2a), 3L))))
  unlink(path)
})

testthat::test_that("multigenerational pedigree BED has zero exact mismatches", {
  backend <- .bed_backend()
  fixture <- .bed_pedigree_fixture(backend)
  path <- tempfile(fileext = ".bed")
  manifest <- .write_one_bed(backend, path, fixture$packed, "chrP", 3L)
  decoded <- gsim:::.gsim_packed_bed_read_all(
    backend, path, length(fixture$packed$sample_ids),
    length(fixture$packed$variant_ids), fixture$packed$sample_ids,
    fixture$packed$variant_ids
  )
  expected <- matrix(as.integer(fixture$raw$genotypes),
                     nrow(fixture$raw$genotypes),
                     dimnames = dimnames(fixture$raw$genotypes))
  testthat::expect_identical(decoded, expected)
  testthat::expect_identical(decoded,
    matrix(as.integer(gsim:::.gsim_packed_decode_genotypes(
      fixture$packed$h1, fixture$packed$h2
    )), nrow(decoded), dimnames = dimnames(decoded)))
  testthat::expect_identical(sum(decoded == -9L), 0L)
  testthat::expect_identical(readBin(path, "raw", n = file.info(path)$size),
                             .bed_expected(fixture$raw$genotypes))
  testthat::expect_identical(manifest$chromosomes$marker_count,
                             ncol(fixture$raw$genotypes))
  unlink(path)
})

testthat::test_that("three chromosomes append exactly in declared order", {
  backend <- .bed_backend()
  labels <- c("chrZ", "01", "CHR1")
  marker_counts <- c(1L, 5L, 3L)

  write_order <- function(order, path, capacity) {
    sink <- gsim:::.gsim_bed_sink_create(
      backend, path, paste0("syn", seq_len(7L)),
      buffer_variants = capacity,
      provenance = list(seed = 410, chromosome_identity = "exact UTF-8 label")
    )
    raw <- list()
    for (label in order) {
      value <- .bed_founder_fixture(
        backend, label, marker_counts[[match(label, labels)]], n = 7L, seed = 410
      )
      gsim:::.gsim_bed_sink_append(
        sink, label, value$packed$h1, value$packed$h2,
        value$packed$variant_ids
      )
      raw[[label]] <- value$raw$genotypes
      rm(value)
      gc(FALSE)
    }
    list(manifest = gsim:::.gsim_bed_sink_finalize(sink),
         expected = do.call(cbind, raw))
  }

  forward1_path <- tempfile(fileext = ".bed")
  forward4_path <- tempfile(fileext = ".bed")
  reverse_path <- tempfile(fileext = ".bed")
  forward1 <- write_order(labels, forward1_path, 1L)
  forward4 <- write_order(labels, forward4_path, 4L)
  reverse <- write_order(rev(labels), reverse_path, 2L)
  expected_size <- 3 + sum(marker_counts) * ceiling(7 / 4)
  testthat::expect_identical(file.info(forward1_path)$size, expected_size)
  testthat::expect_identical(forward1$manifest$bytes_written, expected_size)
  testthat::expect_identical(forward1$manifest$chromosome_order, labels)
  testthat::expect_identical(forward1$manifest$chromosomes$start_variant,
                             c(1L, 2L, 7L))
  testthat::expect_identical(readBin(forward1_path, "raw", n = 999L),
                             .bed_expected(forward1$expected))
  testthat::expect_identical(readBin(forward1_path, "raw", n = 999L),
                             readBin(forward4_path, "raw", n = 999L))
  testthat::expect_identical(readBin(reverse_path, "raw", n = 999L),
                             .bed_expected(reverse$expected))
  testthat::expect_identical(
    gsim:::.gsim_packed_bed_read_all(
      backend, reverse_path, 7L, sum(marker_counts),
      reverse$manifest$sample_ids, reverse$manifest$variant_ids
    ),
    matrix(as.integer(reverse$expected), 7L,
           dimnames = list(reverse$manifest$sample_ids,
                           reverse$manifest$variant_ids))
  )
  testthat::expect_identical(
    sum(readBin(forward1_path, "raw", n = 999L)[1:3] ==
          as.raw(c(0x6c, 0x1b, 0x01))), 3L
  )
  unlink(c(forward1_path, forward4_path, reverse_path))
})

testthat::test_that("R sink validation and transactional filesystem behavior are explicit", {
  backend <- .bed_backend()
  directory <- tempfile("gsim_bed_sink_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  path <- file.path(directory, "output.bed")
  samples <- c("a", "b")
  variants <- "v"
  h1 <- matrix(as.raw(c(0, 1)), 2L, 1L,
               dimnames = list(samples, variants))
  h2 <- matrix(as.raw(c(1, 1)), 2L, 1L, dimnames = dimnames(h1))
  p1 <- gsim:::.gsim_packed_pack(backend, h1)
  p2 <- gsim:::.gsim_packed_pack(backend, h2)

  sink <- gsim:::.gsim_bed_sink_create(backend, path, samples,
                                       buffer_variants = 1L)
  testthat::expect_false(file.exists(path))
  testthat::expect_length(list.files(directory, pattern = "gsim\\.tmp"), 1L)
  gsim:::.gsim_bed_sink_cancel(sink)
  testthat::expect_false(file.exists(path))
  testthat::expect_length(list.files(directory, pattern = "gsim\\.tmp"), 0L)
  testthat::expect_error(gsim:::.gsim_bed_sink_finalize(sink), "cancelled")

  writeBin(as.raw(0x55), path)
  testthat::expect_error(
    gsim:::.gsim_bed_sink_create(backend, path, samples),
    "overwrite is disabled"
  )
  testthat::expect_identical(readBin(path, "raw", n = 9L), as.raw(0x55))
  sink <- gsim:::.gsim_bed_sink_create(backend, path, samples,
                                       overwrite = TRUE)
  gsim:::.gsim_bed_sink_append(sink, "x", p1, p2, variants)
  testthat::expect_identical(readBin(path, "raw", n = 9L), as.raw(0x55))
  manifest <- gsim:::.gsim_bed_sink_finalize(sink)
  testthat::expect_identical(
    readBin(path, "raw", n = 9L),
    .bed_expected(matrix(as.integer(h1) + as.integer(h2), 2L, 1L))
  )
  testthat::expect_error(gsim:::.gsim_bed_sink_finalize(sink),
                         "already been finalized")
  testthat::expect_error(gsim:::.gsim_bed_sink_append(
    sink, "y", p1, p2, variants), "Cannot append")
  testthat::expect_identical(manifest$state, "finalized")

  reordered <- h1[2:1, , drop = FALSE]
  sink2 <- gsim:::.gsim_bed_sink_create(
    backend, file.path(directory, "order.bed"), samples
  )
  testthat::expect_error(gsim:::.gsim_bed_sink_append(
    sink2, "x", gsim:::.gsim_packed_pack(backend, reordered), p2, variants
  ), "sample order")
  gsim:::.gsim_bed_sink_cancel(sink2)
  testthat::expect_error(
    gsim:::.gsim_bed_sink_create(
      backend, file.path(directory, "absent", "bad.bed"), samples
    ), "cannot find|mustWork|parent"
  )
})

testthat::test_that("writer memory is chromosome-local and throughput is linear-scale", {
  backend <- .bed_backend()
  individuals <- 513L
  markers <- 2048L
  h1 <- matrix(as.raw(rep(c(0L, 1L, 1L, 0L),
                          length.out = individuals * markers)),
               individuals, markers)
  h2 <- matrix(as.raw(rep(c(1L, 0L, 1L, 0L),
                          length.out = individuals * markers)),
               individuals, markers)
  dimnames(h1) <- dimnames(h2) <- list(
    paste0("i", seq_len(individuals)), paste0("v", seq_len(markers))
  )
  p1 <- gsim:::.gsim_packed_pack(backend, h1)
  p2 <- gsim:::.gsim_packed_pack(backend, h2)
  path <- tempfile(fileext = ".bed")
  sink <- gsim:::.gsim_bed_sink_create(
    backend, path, rownames(h1), buffer_variants = 64L
  )
  elapsed <- system.time(gsim:::.gsim_bed_sink_append(
    sink, "bounded", p1, p2, colnames(h1)
  ))[["elapsed"]]
  manifest <- gsim:::.gsim_bed_sink_finalize(sink)
  packed_bytes <- unname(gsim:::.gsim_packed_info(p1)[[4L]] +
                           gsim:::.gsim_packed_info(p2)[[4L]])
  expected_packed <- 2 * 8 * markers * ceiling(individuals / 64)
  expected_bed <- 3 + markers * ceiling(individuals / 4)
  expected_buffer <- 64 * ceiling(individuals / 4)
  testthat::expect_identical(packed_bytes, expected_packed)
  testthat::expect_identical(manifest$bytes_written, expected_bed)
  testthat::expect_identical(file.info(path)$size, expected_bed)
  testthat::expect_identical(manifest$conversion_buffer_bytes,
                             expected_buffer)
  testthat::expect_true(is.finite(elapsed) && elapsed < 30)
  testthat::expect_false("genotypes" %in% names(manifest))
  unlink(path)
})
