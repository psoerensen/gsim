.public_pedigree <- function() {
  tab <- data.frame(
    animal = c("F1", "F2", "M1", "M2", "L", "FS1", "FS2", "PHS", "MHS", "G"),
    sire = c(rep(NA_character_, 5), "F1", "F1", "F1", "F2", "FS1"),
    dam = c(rep(NA_character_, 5), "M1", "M1", "M2", "M1", "L"),
    sex = "U", generation = c(rep(1L, 5), rep(2L, 4), 3L),
    cohort = c(rep(1L, 5), rep(2L, 4), 3L), phenotyped = TRUE,
    stringsAsFactors = FALSE
  )
  canonical <- tab$animal
  external <- rev(canonical)
  structure(list(
    pedigree = tab[match(external, canonical), , drop = FALSE],
    canonical_order = canonical, external_order = external,
    mapping = data.frame(animal = canonical,
                         canonical_index = seq_along(canonical),
                         external_index = match(canonical, external)),
    settings = list(), diagnostics = list(), checksums = list()
  ), class = "gsim_pedigree")
}

.public_vcf_write <- function(path, records, samples) {
  header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", samples), collapse = "\t")
  writeBin(charToRaw(paste(c("##fileformat=VCFv4.2", header, records, ""),
                           collapse = "\n")), path)
}

testthat::test_that("public VCF to founder, pedigree, HAP, and BED is exact", {
  root <- tempfile("public-vcf-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  samples <- paste0("d", 1:4)
  records <- c(
    "B\t10\tb1\tA\tG\t.\t.\t.\tGT\t0|1\t0|1\t1|0\t1|0",
    "B\t20\tb2\tC\tT\t.\t.\t.\tDP:GT\t2:1|0\t2:0|1\t2:1|1\t2:0|0",
    "B\t30\tb3\tG\tA\t.\t.\t.\tGT\t0|0\t1|1\t0|1\t1|0",
    "B\t40\tb4\tT\tC\t.\t.\t.\tGT\t1|1\t0|0\t1|0\t0|1",
    "A\t10\ta1\tA\tC\t.\t.\t.\tGT\t0|1\t1|0\t0|0\t1|1",
    "A\t20\ta2\tC\tG\t.\t.\t.\tGT\t1|0\t0|1\t1|1\t0|0",
    "A\t30\ta3\tG\tT\t.\t.\t.\tGT\t0|0\t1|1\t1|0\t0|1",
    "A\t40\ta4\tT\tA\t.\t.\t.\tGT\t1|1\t0|0\t0|1\t1|0"
  )
  vcf <- file.path(root, "reference.vcf")
  .public_vcf_write(vcf, records, samples)
  map <- data.frame(
    chromosome = rep(c("B", "A"), each = 4),
    variant_id = c(paste0("b", 1:4), paste0("a", 1:4)),
    genetic_position_cm = rep(c(0, 50, 100, 150), 2),
    stringsAsFactors = FALSE
  )
  set.seed(12345)
  global_before <- .Random.seed
  reference <- gsim_import_vcf(vcf, map, file.path(root, "reference"))
  testthat::expect_s3_class(reference, "gsim_reference")
  reopened <- gsim_reference(file.path(root, "reference"))
  testthat::expect_identical(reopened$sample_ids, samples)
  populations <- stats::setNames(c("P1", "P1", "P2", "P2"), samples)
  mutation <- stats::setNames(c(1e9, 1e9, 0, 1e9, 1e9, 1e9, 0, 1e9),
                              reference$variant_ids)
  pedigree <- .public_pedigree()
  common <- list(
    reference = reference, pedigree = pedigree, populations = populations,
    ancestry_weights = c(P1 = 0.4, P2 = 0.6), mutation_age = mutation,
    N = c(P1 = 2, P2 = 2), Ne = c(P1 = 4, P2 = 6),
    rho = c(P1 = 0.02, P2 = 0.03), seed = 717
  )
  hap <- do.call(gsim_simulate, c(common, list(
    output = file.path(root, "simulation-hap"), format = "hap")))
  bed <- do.call(gsim_simulate, c(common, list(
    output = file.path(root, "simulation-bed"), format = "bed")))
  testthat::expect_identical(.Random.seed, global_before)
  testthat::expect_identical(hap$sample_ids, pedigree$canonical_order)
  testthat::expect_identical(hap$variant_ids, reference$variant_ids)
  testthat::expect_identical(bed$variant_ids, reference$variant_ids)

  gbits <- gsim:::.gsim_packed_backend()
  gmat <- gsim:::.gsim_metadata_backend()
  output <- gsim:::.gsim_hap_dataset_open(
    gbits, gmat, file.path(root, "simulation-hap"))
  reference_reader <- gsim:::.gsim_hap_dataset_open(
    gbits, gmat, file.path(root, "reference"))
  expected_genotypes <- list()
  for (chromosome in c("B", "A")) {
    rows <- which(reference_reader$variants$chromosome == chromosome)
    input <- gsim:::.gsim_hap_dataset_load_chromosome(reference_reader, chromosome)
    raw_h1 <- gsim:::.gsim_packed_unpack(input$h1)
    raw_h2 <- gsim:::.gsim_packed_unpack(input$h2)
    founder <- gsim:::.gsim_hapnest_founders(
      raw_h1, raw_h2, unname(populations), c(P1 = 0.4, P2 = 0.6),
      c(P1 = 2, P2 = 2), c(P1 = 4, P2 = 6), c(P1 = 0.02, P2 = 0.03),
      reference_reader$variants$genetic_position_cm[rows],
      unname(mutation[reference_reader$variants$variant_id[rows]]),
      5L, 717, rep.int(chromosome, length(rows)), return_segments = FALSE)
    rownames(founder$h1) <- rownames(founder$h2) <- pedigree$canonical_order[1:5]
    expected <- gsim:::.gsim_pedigree_genotypes(
      pedigree, list(h1 = founder$h1, h2 = founder$h2),
      rep.int(chromosome, length(rows)),
      reference_reader$variants$genetic_position_cm[rows] / 100, 717,
      return_crossovers = FALSE)
    actual <- gsim:::.gsim_hap_dataset_load_chromosome(output, chromosome)
    actual_h1 <- gsim:::.gsim_packed_unpack(actual$h1)
    actual_h2 <- gsim:::.gsim_packed_unpack(actual$h2)
    testthat::expect_identical(actual_h1, expected$h1)
    testthat::expect_identical(actual_h2, expected$h2)
    expected_genotypes[[chromosome]] <- expected$genotypes
    for (i in which(!expected$pedigree_alignment$founder)) {
      sire <- match(expected$pedigree_alignment$sire[[i]], expected$sample_ids)
      dam <- match(expected$pedigree_alignment$dam[[i]], expected$sample_ids)
      testthat::expect_true(all(actual_h1[i, ] == actual_h1[sire, ] |
                                  actual_h1[i, ] == actual_h2[sire, ]))
      testthat::expect_true(all(actual_h2[i, ] == actual_h1[dam, ] |
                                  actual_h2[i, ] == actual_h2[dam, ]))
    }
    gsim:::.gsim_packed_close(input$h1); gsim:::.gsim_packed_close(input$h2)
    gsim:::.gsim_packed_close(actual$h1); gsim:::.gsim_packed_close(actual$h2)
  }
  decoded <- gsim:::.gsim_packed_bed_read_all(
    gbits, bed$paths[["bed"]], length(pedigree$canonical_order),
    length(reference$variant_ids))
  expected_all <- do.call(cbind, expected_genotypes[c("B", "A")])
  testthat::expect_identical(as.vector(decoded), as.integer(expected_all))
  testthat::expect_false(anyNA(decoded))

  reverse_vcf <- file.path(root, "reference-reverse.vcf")
  .public_vcf_write(reverse_vcf, c(records[5:8], records[1:4]), samples)
  reverse_reference <- gsim_import_vcf(
    reverse_vcf, map, file.path(root, "reference-reverse"))
  reverse_common <- common
  reverse_common$reference <- reverse_reference
  reverse_hap <- do.call(gsim_simulate, c(reverse_common, list(
    output = file.path(root, "simulation-reverse"), format = "hap")))
  reverse_reader <- gsim:::.gsim_hap_dataset_open(
    gbits, gmat, file.path(root, "simulation-reverse"))
  for (chromosome in c("B", "A")) {
    forward_phase <- gsim:::.gsim_hap_dataset_load_chromosome(output, chromosome)
    reverse_phase <- gsim:::.gsim_hap_dataset_load_chromosome(reverse_reader,
                                                               chromosome)
    testthat::expect_identical(gsim:::.gsim_packed_unpack(forward_phase$h1),
                               gsim:::.gsim_packed_unpack(reverse_phase$h1))
    testthat::expect_identical(gsim:::.gsim_packed_unpack(forward_phase$h2),
                               gsim:::.gsim_packed_unpack(reverse_phase$h2))
    gsim:::.gsim_packed_close(forward_phase$h1)
    gsim:::.gsim_packed_close(forward_phase$h2)
    gsim:::.gsim_packed_close(reverse_phase$h1)
    gsim:::.gsim_packed_close(reverse_phase$h2)
  }
  gsim:::.gsim_hap_dataset_close(reverse_reader)

  single_vcf <- file.path(root, "reference-single.vcf")
  .public_vcf_write(single_vcf, records[1:4], samples)
  single_reference <- gsim_import_vcf(
    single_vcf, map[map$chromosome == "B", ], file.path(root, "reference-single"))
  single_common <- common
  single_common$reference <- single_reference
  single_common$mutation_age <- mutation[paste0("b", 1:4)]
  single_hap <- do.call(gsim_simulate, c(single_common, list(
    output = file.path(root, "simulation-single"), format = "hap")))
  single_reader <- gsim:::.gsim_hap_dataset_open(
    gbits, gmat, file.path(root, "simulation-single"))
  collection_b <- gsim:::.gsim_hap_dataset_load_chromosome(output, "B")
  standalone_b <- gsim:::.gsim_hap_dataset_load_chromosome(single_reader, "B")
  testthat::expect_identical(gsim:::.gsim_packed_unpack(collection_b$h1),
                             gsim:::.gsim_packed_unpack(standalone_b$h1))
  testthat::expect_identical(gsim:::.gsim_packed_unpack(collection_b$h2),
                             gsim:::.gsim_packed_unpack(standalone_b$h2))
  gsim:::.gsim_packed_close(collection_b$h1)
  gsim:::.gsim_packed_close(collection_b$h2)
  gsim:::.gsim_packed_close(standalone_b$h1)
  gsim:::.gsim_packed_close(standalone_b$h2)
  gsim:::.gsim_hap_dataset_close(single_reader)
  gsim:::.gsim_hap_dataset_close(output)
  gsim:::.gsim_hap_dataset_close(reference_reader)

  repeat_hap <- do.call(gsim_simulate, c(common, list(
    output = file.path(root, "simulation-repeat"), format = "hap")))
  testthat::expect_identical(readBin(hap$paths[["hap"]], "raw", 100000),
                             readBin(repeat_hap$paths[["hap"]], "raw", 100000))
  changed <- common
  changed$seed <- 718
  changed_hap <- do.call(gsim_simulate, c(changed, list(
    output = file.path(root, "simulation-changed"), format = "hap")))
  testthat::expect_false(identical(readBin(hap$paths[["hap"]], "raw", 100000),
                                   readBin(changed_hap$paths[["hap"]], "raw", 100000)))
})
