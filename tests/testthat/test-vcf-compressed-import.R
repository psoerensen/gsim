.vcfgz_text <- function(lines) paste0(paste(lines, collapse = "\n"), "\n")

.vcfgz_write_gzip <- function(path, text) {
  connection <- gzfile(path, "wb")
  on.exit(close(connection), add = TRUE)
  writeBin(charToRaw(text), connection)
  invisible(path)
}

.vcfgz_bytes <- function(path) {
  readBin(path, "raw", n = file.info(path)$size)
}

.vcfgz_load <- function(prefix, chromosome) {
  packed <- gsim:::.gsim_packed_backend()
  metadata <- gsim:::.gsim_metadata_backend()
  reader <- gsim:::.gsim_hap_dataset_open(packed, metadata, prefix)
  on.exit(gsim:::.gsim_hap_dataset_close(reader), add = TRUE)
  loaded <- gsim:::.gsim_hap_dataset_load_chromosome(reader, chromosome)
  on.exit({
    gsim:::.gsim_packed_close(loaded$h1)
    gsim:::.gsim_packed_close(loaded$h2)
  }, add = TRUE)
  list(
    h1 = gsim:::.gsim_packed_unpack(loaded$h1),
    h2 = gsim:::.gsim_packed_unpack(loaded$h2),
    h1_words = gsim:::.gsim_packed_word(loaded$h1, 1L, 1L),
    h2_words = gsim:::.gsim_packed_word(loaded$h2, 1L, 1L)
  )
}

testthat::test_that("plain, gzip, and concatenated gzip imports are exactly equal", {
  root <- tempfile("vcfgz-parity-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  header <- c(
    "##fileformat=VCFv4.2",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", "s1", "s2", "s3"), collapse = "\t")
  )
  records <- c(
    "x\t1\tx1\tA\tG\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1",
    "z\t5\tz5\tA\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1",
    "z\t10\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0|1:2\t1|1:3\t1|0:4",
    "z\t15\tz15\tC\tT\t.\tPASS\t.\tDP:GT\t2:1|0\t3:1|1\t4:0|1",
    "z\t20\tz20\tG\tA\t.\tPASS\t.\tGT\t1|1\t0|0\t1|0",
    "z\t25\tz25\tT\tC\t.\tPASS\t.\tGT\t0|0\t0|1\t1|1"
  )
  plain <- file.path(root, "reference.vcf")
  writeBin(charToRaw(.vcfgz_text(c(header, records))), plain)
  gzip <- file.path(root, "reference.vcf.gz")
  .vcfgz_write_gzip(gzip, .vcfgz_text(c(header, records)))
  first_member <- file.path(root, "first.gz")
  second_member <- file.path(root, "second.gz")
  .vcfgz_write_gzip(first_member, .vcfgz_text(c(header, records[1:2])))
  .vcfgz_write_gzip(second_member, .vcfgz_text(records[3:6]))
  concatenated <- file.path(root, "concatenated.data")
  writeBin(c(.vcfgz_bytes(first_member), .vcfgz_bytes(second_member)), concatenated)

  map <- data.frame(
    chromosome = c("z", "z"), base_pair_position = c(10, 20),
    genetic_position_cm = c(0, 1), stringsAsFactors = FALSE
  )
  sources <- c(plain = plain, gzip = gzip, concatenated = concatenated)
  results <- lapply(names(sources), function(name) {
    gsim_import_vcf(
      sources[[name]], map, file.path(root, paste0("out-", name)),
      samples = c("s3", "s1"), chromosome = "z", region = c(10, 20),
      unsupported = "skip"
    )
  })
  names(results) <- names(sources)
  expected_h1 <- matrix(as.raw(c(1, 0, 0, 1, 1, 1)), 2L, 3L)
  expected_h2 <- matrix(as.raw(c(0, 1, 1, 0, 0, 1)), 2L, 3L)
  loaded <- lapply(results, function(value) .vcfgz_load(value$prefix, "z"))
  for (value in loaded) {
    testthat::expect_identical(unname(value$h1), expected_h1)
    testthat::expect_identical(unname(value$h2), expected_h2)
  }
  for (extension in c("hap", "bim", "fam")) {
    bytes <- lapply(results, function(value) .vcfgz_bytes(value$paths[[extension]]))
    testthat::expect_identical(bytes[[1L]], bytes[[2L]])
    testthat::expect_identical(bytes[[1L]], bytes[[3L]])
  }
  testthat::expect_identical(results$plain$sample_ids, c("s3", "s1"))
  testthat::expect_identical(results$plain$variant_ids,
                             c("z:10:A:G", "z15", "z20"))
  testthat::expect_identical(
    readLines(results$plain$paths[["bim"]], warn = FALSE),
    c("z\tz:10:A:G\t0\t10\tG\tA", "z\tz15\t0.5\t15\tT\tC",
      "z\tz20\t1\t20\tA\tG")
  )
  testthat::expect_identical(
    readLines(results$plain$paths[["fam"]], warn = FALSE),
    c("reference\ts3\t0\t0\t0\t-9", "reference\ts1\t0\t0\t0\t-9")
  )
  testthat::expect_identical(results$plain$import$input_type, "plain VCF")
  testthat::expect_identical(results$gzip$import$input_type, "gzip VCF")
  testthat::expect_identical(results$concatenated$import$input_type, "gzip VCF")
  report <- results$gzip$import
  testthat::expect_equal(report$total_records_scanned, 6)
  testthat::expect_equal(report$retained_variants, 3)
  testthat::expect_equal(report$outside_selected_chromosome, 1)
  testthat::expect_equal(report$outside_selected_region, 2)
  testthat::expect_false(report$dense_haplotype_matrix_allocated)
  testthat::expect_false(report$dense_genotype_matrix_allocated)
  testthat::expect_false(report$full_uncompressed_temporary_vcf_created)
  testthat::expect_gt(report$maximum_parsing_buffer_bytes, 65536)
})

testthat::test_that("unsupported records are classified and malformed compression fails", {
  root <- tempfile("vcfgz-filter-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  header <- c(
    "##fileformat=VCFv4.2",
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", "s1", "s2"), collapse = "\t")
  )
  records <- c(
    "1\t1\tv1\tA\tG\t.\t.\t.\tGT\t0|1\tnot-decoded",
    "1\t2\tv2\tA\tAT\t.\t.\t.\tGT\t0|1\tnot-decoded",
    "1\t3\tv3\tA\tG,T\t.\t.\t.\tGT\t0|1\tnot-decoded",
    "1\t4\tv4\tA\t<DEL>\t.\t.\t.\tGT\t0|1\tnot-decoded",
    "1\t5\tv5\tA\tg\t.\t.\t.\tGT\t0|1\tnot-decoded",
    "1\t6\tv6\tA\tG\t.\t.\t.\tGT\t.|1\tnot-decoded",
    "1\t7\tv7\tA\tG\t.\t.\t.\tGT\t0/1\tnot-decoded",
    "1\t8\tv8\tA\tG\t.\t.\t.\tGT\t1\tnot-decoded"
  )
  plain <- file.path(root, "mixed.vcf")
  writeBin(charToRaw(.vcfgz_text(c(header, records))), plain)
  map <- data.frame(chromosome = c("1", "1"),
                    base_pair_position = c(1, 8),
                    genetic_position_cm = c(0, 0.7))
  result <- gsim_import_vcf(
    plain, map, file.path(root, "skip"), samples = "s1",
    chromosome = "1", unsupported = "skip"
  )
  report <- result$import
  testthat::expect_equal(report$retained_variants, 1)
  testthat::expect_equal(unname(unlist(report[c(
    "indels", "multiallelic_records", "symbolic_or_breakend_alleles",
    "other_unsupported_alleles", "missing_gt", "unphased_gt",
    "non_diploid_gt")])), rep(1, 7))
  testthat::expect_error(gsim_import_vcf(
    plain, map, file.path(root, "error"), samples = "s1",
    chromosome = "1", unsupported = "error"), "1:2.*v2.*indel")
  testthat::expect_error(gsim_import_vcf(
    plain, map, file.path(root, "unknown"), samples = "absent",
    chromosome = "1"), "absent")
  testthat::expect_error(gsim_import_vcf(
    plain, map, file.path(root, "duplicate-samples"), samples = c("s1", "s1"),
    chromosome = "1"), "unique")

  malformed <- file.path(root, "malformed.vcf")
  writeBin(charToRaw(.vcfgz_text(c(header,
    "1\t1\tv\tA\tG\t.\t.\t.\tGT\t0||1\t0|0"))), malformed)
  testthat::expect_error(gsim_import_vcf(
    malformed, map, file.path(root, "malformed"), samples = "s1",
    chromosome = "1", unsupported = "skip"), "malformed GT")

  valid_gzip <- file.path(root, "valid.gz")
  .vcfgz_write_gzip(valid_gzip, .vcfgz_text(c(header, records[[1L]])))
  compressed <- .vcfgz_bytes(valid_gzip)
  truncated <- file.path(root, "truncated.gz")
  writeBin(head(compressed, -4L), truncated)
  testthat::expect_error(gsim_import_vcf(
    truncated, map, file.path(root, "truncated"), samples = "s1",
    chromosome = "1", unsupported = "skip"), "truncated|corrupt|checksum")

  outside_map <- data.frame(chromosome = c("1", "1"),
                            base_pair_position = c(2, 8),
                            genetic_position_cm = c(0, 1))
  testthat::expect_error(gsim_import_vcf(
    plain, outside_map, file.path(root, "outside-map"), samples = "s1",
    chromosome = "1", unsupported = "skip"), "outside.*map range")
})
