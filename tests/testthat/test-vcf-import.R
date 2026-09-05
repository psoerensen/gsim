.vcf_backends <- function() {
  list(
    gbits = gsim:::.gsim_packed_backend(),
    gmat = gsim:::.gsim_metadata_backend()
  )
}

.vcf_write <- function(path, records, samples = c("s1", "s2")) {
  header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", samples), collapse = "\t")
  writeBin(charToRaw(paste(c("##fileformat=VCFv4.2", header, records, ""),
                           collapse = "\n")), path)
}

testthat::test_that("strict VCF import preserves exact phase and metadata", {
  backend <- .vcf_backends()
  root <- tempfile("vcf-import-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  vcf <- file.path(root, "reference.vcf")
  .vcf_write(vcf, c(
    "z\t10\t.\tA\tG\t.\tPASS\t.\tDP:GT\t5:0|1\t7:1|0",
    "z\t12\tv2\tC\tT\t.\tPASS\t.\tGT\t1|1\t0|0",
    "2\t20\tv3\tG\tA\t.\tPASS\t.\tGT:DP\t0|0:3\t1|1:4"
  ))
  map <- data.frame(
    chromosome = c("2", "z", "z"),
    base_pair_position = c(20, 12, 10),
    genetic_position_cm = c(0, 0.25, 0), stringsAsFactors = FALSE
  )
  samples <- data.frame(
    individual_id = c("s2", "s1"), family_id = c("f2", "f1"),
    sex = c(2L, 1L), stringsAsFactors = FALSE
  )
  set.seed(991)
  before <- .Random.seed
  first <- gsim:::.gsim_import_vcf_internal(
    backend$gbits, backend$gmat, vcf, map, file.path(root, "reference"),
    samples
  )
  testthat::expect_identical(.Random.seed, before)
  testthat::expect_identical(first$sample_ids, c("s1", "s2"))
  testthat::expect_identical(first$variant_ids, c("z:10:A:G", "v2", "v3"))
  testthat::expect_identical(first$import$generated_variant_ids, "z:10:A:G")
  reader <- gsim:::.gsim_hap_dataset_open(
    backend$gbits, backend$gmat, file.path(root, "reference"))
  z <- gsim:::.gsim_hap_dataset_load_chromosome(reader, "z")
  two <- gsim:::.gsim_hap_dataset_load_chromosome(reader, "2")
  testthat::expect_identical(
    unname(gsim:::.gsim_packed_unpack(z$h1)),
    matrix(as.raw(c(0, 1, 1, 0)), 2L, 2L))
  testthat::expect_identical(
    unname(gsim:::.gsim_packed_unpack(z$h2)),
    matrix(as.raw(c(1, 0, 1, 0)), 2L, 2L))
  testthat::expect_identical(
    unname(gsim:::.gsim_packed_unpack(two$h1)), matrix(as.raw(c(0, 1)), 2L, 1L))
  testthat::expect_identical(
    unname(gsim:::.gsim_packed_unpack(two$h2)), matrix(as.raw(c(0, 1)), 2L, 1L))
  testthat::expect_identical(
    as.integer(gsim:::.gsim_packed_word(z$h1, 1L, 1L)),
    c(2L, rep.int(0L, 7L)))
  testthat::expect_identical(
    as.integer(gsim:::.gsim_packed_word(z$h2, 1L, 1L)),
    c(1L, rep.int(0L, 7L)))
  testthat::expect_identical(
    readLines(first$paths[["bim"]], warn = FALSE),
    c("z\tz:10:A:G\t0\t10\tG\tA", "z\tv2\t0.25\t12\tT\tC",
      "2\tv3\t0\t20\tA\tG"))
  testthat::expect_identical(
    readLines(first$paths[["fam"]], warn = FALSE),
    c("f1\ts1\t0\t0\t1\t-9", "f2\ts2\t0\t0\t2\t-9"))
  gsim:::.gsim_packed_close(z$h1); gsim:::.gsim_packed_close(z$h2)
  gsim:::.gsim_packed_close(two$h1); gsim:::.gsim_packed_close(two$h2)
  gsim:::.gsim_hap_dataset_close(reader)

  second <- gsim:::.gsim_import_vcf_internal(
    backend$gbits, backend$gmat, vcf, map, file.path(root, "repeat"), samples
  )
  testthat::expect_identical(readBin(first$paths[["hap"]], "raw", 10000),
                             readBin(second$paths[["hap"]], "raw", 10000))
  testthat::expect_equal(first$observed_hap_bytes,
                         64 + 48 * 2 + 2 * 3 * 8)
})

testthat::test_that("packed import handles word-boundary sample counts", {
  backend <- .vcf_backends()
  root <- tempfile("vcf-boundary-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (n in c(1L, 63L, 64L, 65L)) {
    ids <- paste0("s", seq_len(n))
    gt <- rep(c("0|0", "0|1", "1|0", "1|1"), length.out = n)
    vcf <- file.path(root, paste0("n", n, ".vcf"))
    .vcf_write(vcf, paste(c("1", "1", "v", "A", "C", ".", ".", ".",
                            "GT", gt), collapse = "\t"), ids)
    manifest <- gsim:::.gsim_import_vcf_internal(
      backend$gbits, backend$gmat, vcf,
      data.frame(chromosome = "1", variant_id = "v",
                 genetic_position_cm = 0), file.path(root, paste0("out", n)))
    reader <- gsim:::.gsim_hap_dataset_open(
      backend$gbits, backend$gmat, sub("\\.hap$", "", manifest$paths[["hap"]]))
    loaded <- gsim:::.gsim_hap_dataset_load_chromosome(reader, "1")
    expected_h1 <- as.raw(as.integer(substr(gt, 1, 1)))
    expected_h2 <- as.raw(as.integer(substr(gt, 3, 3)))
    testthat::expect_identical(as.vector(gsim:::.gsim_packed_unpack(loaded$h1)),
                               expected_h1)
    testthat::expect_identical(as.vector(gsim:::.gsim_packed_unpack(loaded$h2)),
                               expected_h2)
    testthat::expect_equal(manifest$packed_data_bytes, 2 * ceiling(n / 64) * 8)
    gsim:::.gsim_packed_close(loaded$h1); gsim:::.gsim_packed_close(loaded$h2)
    gsim:::.gsim_hap_dataset_close(reader)
  }
})

testthat::test_that("strict VCF parser rejects every excluded call class", {
  backend <- .vcf_backends()
  root <- tempfile("vcf-invalid-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  bad <- list(
    unphased = "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0/1",
    missing = "1\t1\tv\tA\tG\t.\t.\t.\tGT\t.|1",
    haploid = "1\t1\tv\tA\tG\t.\t.\t.\tGT\t1",
    polyploid = "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0|1|0",
    multiallelic = "1\t1\tv\tA\tG,T\t.\t.\t.\tGT\t0|1",
    symbolic = "1\t1\tv\tA\t<DEL>\t.\t.\t.\tGT\t0|1",
    breakend = "1\t1\tv\tA\tA]2:3]\t.\t.\t.\tGT\t0|1",
    indel = "1\t1\tv\tA\tAT\t.\t.\t.\tGT\t0|1",
    lowercase = "1\t1\tv\ta\tG\t.\t.\t.\tGT\t0|1",
    identical_alleles = "1\t1\tv\tA\tA\t.\t.\t.\tGT\t0|1",
    no_gt = "1\t1\tv\tA\tG\t.\t.\t.\tDP\t2",
    duplicate_gt = "1\t1\tv\tA\tG\t.\t.\t.\tGT:GT\t0|1:0|1",
    bad_arity = "1\t1\tv\tA\tG\t.\t.\t.\tGT:DP\t0|1",
    bad_pos = "1\tx\tv\tA\tG\t.\t.\t.\tGT\t0|1",
    empty_id = "1\t1\t\tA\tG\t.\t.\t.\tGT\t0|1"
  )
  for (name in names(bad)) {
    path <- file.path(root, paste0(name, ".vcf"))
    .vcf_write(path, bad[[name]], "s")
    testthat::expect_error(gsim:::.gsim_vcf_reader_open(backend$gmat, path),
                           "VCF", info = name)
  }
  duplicate_samples <- file.path(root, "dup-samples.vcf")
  .vcf_write(duplicate_samples,
             "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0|0\t1|1", c("s", "s"))
  testthat::expect_error(
    gsim:::.gsim_vcf_reader_open(backend$gmat, duplicate_samples), "sample IDs")
  disjoint <- file.path(root, "disjoint.vcf")
  .vcf_write(disjoint, c(
    "1\t1\tv1\tA\tG\t.\t.\t.\tGT\t0|0",
    "2\t1\tv2\tA\tG\t.\t.\t.\tGT\t0|0",
    "1\t2\tv3\tA\tG\t.\t.\t.\tGT\t0|0"), "s")
  testthat::expect_error(gsim:::.gsim_vcf_reader_open(backend$gmat, disjoint),
                         "disjoint")
  duplicate_id <- file.path(root, "dup-id.vcf")
  .vcf_write(duplicate_id, c(
    "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0|0",
    "1\t2\tv\tC\tT\t.\t.\t.\tGT\t0|0"), "s")
  testthat::expect_error(gsim:::.gsim_vcf_reader_open(backend$gmat, duplicate_id),
                         "duplicate final variant")
  wrong_fields <- file.path(root, "wrong-fields.vcf")
  .vcf_write(wrong_fields, "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0|0", c("s1", "s2"))
  testthat::expect_error(gsim:::.gsim_vcf_reader_open(backend$gmat, wrong_fields),
                         "sample-field count")
  missing_fileformat <- file.path(root, "missing-fileformat.vcf")
  writeBin(charToRaw(paste0(
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts\n",
    "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0|0\n")), missing_fileformat)
  testthat::expect_error(
    gsim:::.gsim_vcf_reader_open(backend$gmat, missing_fileformat), "fileformat")
  compressed <- file.path(root, "unsupported.vcf.gz")
  file.copy(duplicate_id, compressed)
  testthat::expect_error(gsim:::.gsim_vcf_reader_open(backend$gmat, compressed),
                         "unsupported")
})

testthat::test_that("map and publication failures leave no partial dataset", {
  backend <- .vcf_backends()
  root <- tempfile("vcf-cleanup-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  vcf <- file.path(root, "x.vcf")
  .vcf_write(vcf, "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0|1", "s")
  prefix <- file.path(root, "bad")
  testthat::expect_error(gsim:::.gsim_import_vcf_internal(
    backend$gbits, backend$gmat, vcf,
    data.frame(chromosome = "1", variant_id = "other",
               genetic_position_cm = 0), prefix), "alignment")
  testthat::expect_false(any(file.exists(paste0(prefix, c(".hap", ".bim", ".fam")))))
  testthat::expect_error(gsim:::.gsim_import_vcf_internal(
    backend$gbits, backend$gmat, vcf,
    data.frame(chromosome = "1", variant_id = "v",
               genetic_position_cm = -1), prefix), "nonnegative")
  testthat::expect_false(any(file.exists(paste0(prefix, c(".hap", ".bim", ".fam")))))
  testthat::expect_error(gsim:::.gsim_import_vcf_internal(
    backend$gbits, backend$gmat, vcf,
    data.frame(chromosome = "1", variant_id = "v",
               genetic_position_cm = 0), prefix,
    data.frame(individual_id = c("s", "extra"))), "exactly once")
  testthat::expect_false(any(file.exists(paste0(prefix, c(".hap", ".bim", ".fam")))))
})
