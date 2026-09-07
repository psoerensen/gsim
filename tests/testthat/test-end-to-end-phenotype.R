.e2e_write_reference <- function(path, donor_ids) {
  chromosomes <- c("2", "1")
  markers_per_chromosome <- 16L
  variant_ids <- unlist(lapply(chromosomes, function(chromosome) {
    paste0("chr", chromosome, "_m", seq_len(markers_per_chromosome))
  }), use.names = FALSE)
  header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", donor_ids), collapse = "\t")
  patterns <- c("0|0", "0|1", "1|0", "1|1")
  records <- unlist(lapply(seq_along(chromosomes), function(block) {
    chromosome <- chromosomes[[block]]
    vapply(seq_len(markers_per_chromosome), function(marker) {
      calls <- patterns[
        ((seq_along(donor_ids) + marker + block - 3L) %% 4L) + 1L
      ]
      paste(c(chromosome, marker * 100L,
              paste0("chr", chromosome, "_m", marker),
              "A", "G", ".", "PASS", ".", "GT", calls), collapse = "\t")
    }, character(1L))
  }), use.names = FALSE)
  writeLines(c("##fileformat=VCFv4.2", header, records), path, useBytes = TRUE)
  list(
    variant_ids = variant_ids,
    map = data.frame(
      chromosome = rep(chromosomes, each = markers_per_chromosome),
      variant_id = variant_ids,
      genetic_position_cm = rep(
        seq(0, 150, length.out = markers_per_chromosome),
        length(chromosomes)
      ),
      stringsAsFactors = FALSE
    )
  )
}

.e2e_read_phases <- function(prefix) {
  reader <- gsim:::.gsim_hap_dataset_open(prefix)
  on.exit(gsim:::.gsim_hap_dataset_close(reader), add = TRUE)
  phases <- setNames(lapply(reader$chromosome, function(chromosome) {
    packed <- gsim:::.gsim_hap_dataset_load_chromosome(reader, chromosome)
    on.exit({
      gsim:::.gsim_packed_close(packed$h1)
      gsim:::.gsim_packed_close(packed$h2)
    }, add = TRUE)
    list(
      h1 = gsim:::.gsim_packed_unpack(packed$h1),
      h2 = gsim:::.gsim_packed_unpack(packed$h2)
    )
  }), reader$chromosome)
  list(phases = phases, samples = reader$samples$individual_id,
       variants = reader$variants)
}

testthat::test_that("packed pedigree BED enters all four Glist phenotype models", {
  testthat::skip_if_not_installed("qgg")
  root <- tempfile("end-to-end-phenotype-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  donor_ids <- paste0("donor", seq_len(8L))
  vcf <- file.path(root, "reference.vcf")
  fixture <- .e2e_write_reference(vcf, donor_ids)
  reference <- gsim_import_vcf(
    vcf, fixture$map, file.path(root, "reference"), unsupported = "error"
  )
  reference <- gsim_reference(reference$prefix)

  pedigree <- gsim_pedigree(
    n_generations = 3L, animals_per_generation = 8L,
    founder_generations = 1L, sires_per_generation = 2L,
    dams_per_generation = 3L, overlapping_generation_probability = 0,
    unknown_sire_probability = 0, unknown_dam_probability = 0,
    new_founder_probability = 0.15, unphenotyped_probability = 0,
    seed = 2026
  )
  testthat::expect_gt(pedigree$diagnostics$full_sib_families, 0)
  testthat::expect_gt(pedigree$diagnostics$paternal_half_sib_sires, 0)
  testthat::expect_gt(pedigree$diagnostics$maternal_half_sib_dams, 0)
  aligned_pedigree <- pedigree$pedigree[
    match(pedigree$canonical_order, pedigree$pedigree$animal), , drop = FALSE
  ]
  founder_ids <- aligned_pedigree$animal[
    is.na(aligned_pedigree$sire) & is.na(aligned_pedigree$dam)
  ]
  populations <- stats::setNames(
    rep(c("P1", "P2"), each = 4L), donor_ids
  )
  mutation_age <- stats::setNames(
    rep(1e9, length(fixture$variant_ids)), fixture$variant_ids
  )
  base <- gsim_simulate_founders(
    reference = reference, founder_ids = founder_ids,
    populations = populations, ancestry_weights = c(P1 = 0.5, P2 = 0.5),
    mutation_age = mutation_age, N = c(P1 = 4, P2 = 4),
    Ne = c(P1 = 100, P2 = 120), rho = c(P1 = 0.02, P2 = 0.03),
    seed = 7001, output = file.path(root, "base"),
    batch_size = 64L, threads = 2L
  )
  pedigree_hap_prefix <- file.path(root, "pedigree-hap")
  pedigree_hap <- gsim_simulate_pedigree(
    base, pedigree, seed = 7002, output = pedigree_hap_prefix,
    format = "hap"
  )
  pedigree_bed <- gsim_simulate_pedigree(
    base, pedigree, seed = 7002, output = file.path(root, "pedigree-bed"),
    format = "bed"
  )

  phased <- .e2e_read_phases(pedigree_hap_prefix)
  expected_dosage <- do.call(cbind, lapply(phased$phases, function(x) {
    matrix(as.integer(x$h1) + as.integer(x$h2), nrow = nrow(x$h1))
  }))
  bed_dosage <- gsim:::.gsim_packed_bed_read_all(
    pedigree_bed$paths[["bed"]], length(pedigree$canonical_order),
    length(fixture$variant_ids)
  )
  testthat::expect_identical(unname(bed_dosage), unname(expected_dosage))
  testthat::expect_identical(phased$samples, pedigree$canonical_order)
  testthat::expect_identical(
    phased$variants$variant_id, fixture$variant_ids
  )

  mendelian_inconsistencies <- 0L
  for (phase in phased$phases) {
    for (i in which(!is.na(aligned_pedigree$sire))) {
      sire <- match(aligned_pedigree$sire[[i]], pedigree$canonical_order)
      dam <- match(aligned_pedigree$dam[[i]], pedigree$canonical_order)
      mendelian_inconsistencies <- mendelian_inconsistencies +
        sum(!(phase$h1[i, ] == phase$h1[sire, ] |
              phase$h1[i, ] == phase$h2[sire, ])) +
        sum(!(phase$h2[i, ] == phase$h1[dam, ] |
              phase$h2[i, ] == phase$h2[dam, ]))
    }
  }
  testthat::expect_equal(mendelian_inconsistencies, 0L)

  Glist <- suppressMessages(qgg::gprep(
    study = "gsim-end-to-end-test",
    bedfiles = pedigree_bed$paths[["bed"]],
    bimfiles = pedigree_bed$paths[["bim"]],
    famfiles = pedigree_bed$paths[["fam"]]
  ))
  marker_ids <- as.character(unlist(Glist$rsids, use.names = FALSE))
  testthat::expect_identical(Glist$ids, pedigree$canonical_order)
  testthat::expect_identical(marker_ids, fixture$variant_ids)
  glist_dosage <- qgg::getG(
    Glist = Glist, chr = 1L, rsids = marker_ids,
    ids = pedigree$canonical_order, impute = TRUE, scale = FALSE
  )
  testthat::expect_identical(
    as.integer(glist_dosage), as.integer(bed_dosage)
  )
  testthat::expect_identical(rownames(glist_dosage), pedigree$canonical_order)
  testthat::expect_identical(colnames(glist_dosage), marker_ids)

  requested <- new.env(parent = emptyenv())
  requested$rsids <- character(0)
  tracking_getG <- function(Glist, rsids, ids, chr = NULL,
                            impute = TRUE, scale = FALSE) {
    requested$rsids <- c(requested$rsids, rsids)
    qgg::getG(
      Glist = Glist, chr = chr, rsids = rsids, ids = ids,
      impute = impute, scale = scale
    )
  }
  common <- list(
    Glist = Glist, architecture = "bayesr", standardize_W = FALSE,
    scale_effects = FALSE, return_genotypes = FALSE,
    getG_fun = tracking_getG
  )
  default <- do.call(gsim, c(common, list(n_causal = 8L, seed = 8001)))
  rng_default <- .Random.seed
  explicit_default <- do.call(gsim, c(common, list(
    n_causal = 8L, seed = 8001, a = 0, b = 0, c = 0,
    ld_score = "unused", annotation_score = "unused"
  )))
  rng_explicit_default <- .Random.seed
  for (field in c("causal_rsids", "component", "B", "B_causal", "G",
                  "E", "Y", "causal_probability", "marker_multipliers")) {
    testthat::expect_identical(default[[field]], explicit_default[[field]])
  }
  testthat::expect_identical(rng_default, rng_explicit_default)
  testthat::expect_identical(rownames(default$Y), pedigree$canonical_order)
  testthat::expect_equal(
    unname(default$G),
    unname(bed_dosage[, match(default$causal_rsids, marker_ids), drop = FALSE] %*% default$B_causal),
    tolerance = 1e-12
  )
  testthat::expect_setequal(unique(requested$rsids), default$causal_rsids)
  testthat::expect_lt(length(unique(requested$rsids)), length(marker_ids))

  maf <- stats::setNames(
    as.numeric(unlist(Glist$maf, use.names = FALSE)), marker_ids
  )
  ld_score <- stats::setNames(
    seq(0.8, 2, length.out = length(marker_ids)), marker_ids
  )
  annotation_score <- stats::setNames(
    seq(0.75, 1.25, length.out = length(marker_ids)), marker_ids
  )
  q <- stats::setNames(rep(c(1, 0, 0, 0), length.out = length(marker_ids)),
                       marker_ids)
  q_only <- do.call(gsim, c(common, list(
    causal_probability = q[rev(names(q))], seed = 8002
  )))
  variance <- do.call(gsim, c(common, list(
    n_causal = 8L, a = -0.4, b = -1,
    ld_score = ld_score[rev(names(ld_score))], seed = 8001
  )))
  combined <- do.call(gsim, c(common, list(
    causal_probability = q, a = -0.4, b = -1, c = 0.5,
    ld_score = ld_score, annotation_score = annotation_score,
    seed = 8002
  )))
  expected_probability <- cbind(
    1 - q, q * 0.8, q * 0.18, q * 0.02
  )
  dimnames(expected_probability) <- dimnames(q_only$marker_probabilities)
  expected_variance <- (maf * (1 - maf))^-0.4 * ld_score^-1
  expected_combined <- expected_variance * annotation_score^0.5

  testthat::expect_equal(q_only$marker_probabilities, expected_probability,
                         tolerance = 1e-15)
  testthat::expect_identical(q_only$causal_probability, q)
  testthat::expect_true(all(q_only$marker_multipliers == 1))
  testthat::expect_identical(default$component, variance$component)
  testthat::expect_identical(default$marker_probabilities,
                             variance$marker_probabilities)
  testthat::expect_equal(variance$marker_multipliers, expected_variance,
                         tolerance = 1e-15)
  testthat::expect_equal(
    variance$B[default$component != 1L, , drop = FALSE],
    default$B[default$component != 1L, , drop = FALSE] *
      sqrt(expected_variance[default$component != 1L]),
    tolerance = 1e-15
  )
  testthat::expect_identical(q_only$component, combined$component)
  testthat::expect_identical(q_only$marker_probabilities,
                             combined$marker_probabilities)
  testthat::expect_identical(combined$causal_probability, q)
  testthat::expect_equal(combined$marker_multipliers, expected_combined,
                         tolerance = 1e-15)
  active <- q_only$component != 1L
  testthat::expect_equal(
    combined$B[active, , drop = FALSE],
    q_only$B[active, , drop = FALSE] * sqrt(expected_combined[active]),
    tolerance = 1e-15
  )
  testthat::expect_equal(
    combined$causal$causal_probability,
    unname(q[combined$causal$rsid]), tolerance = 0
  )
  testthat::expect_equal(
    combined$causal$variance_weight,
    unname(expected_combined[combined$causal$rsid]), tolerance = 1e-15
  )
  testthat::expect_identical(
    combined$settings$marker_multipliers$sources$maf$source, "Glist"
  )
  testthat::expect_identical(
    combined$settings$marker_multipliers$sources$ld_score$source,
    "explicit_input"
  )
  testthat::expect_identical(
    combined$settings$marker_multipliers$sources$annotation_score$source,
    "explicit_input"
  )

  causal_matrix_bytes <- as.numeric(utils::object.size(bed_dosage[, match(default$causal_rsids, marker_ids), drop = FALSE]))
  full_diagnostic_bytes <- as.numeric(utils::object.size(glist_dosage))
  testthat::expect_lt(causal_matrix_bytes, full_diagnostic_bytes)
  testthat::expect_null(default$W_causal)
  testthat::expect_equal(ncol(glist_dosage), 32L)
})
