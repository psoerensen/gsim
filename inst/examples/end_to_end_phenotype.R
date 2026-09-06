# End-to-end packed genotype simulation and Glist phenotype simulation.
# qgg is optional for gsim generally, but is required for this Glist example.
if (!requireNamespace("qgg", quietly = TRUE)) {
  stop("This example requires the suggested qgg package.")
}

local({
  root <- tempfile("gsim-end-to-end-phenotype-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  donor_ids <- paste0("donor", seq_len(8L))
  chromosomes <- c("2", "1")
  markers_per_chromosome <- 16L
  variant_ids <- unlist(lapply(chromosomes, function(chromosome) {
    paste0("chr", chromosome, "_m", seq_len(markers_per_chromosome))
  }), use.names = FALSE)

  header <- paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", donor_ids), collapse = "\t")
  gt <- c("0|0", "0|1", "1|0", "1|1")
  records <- unlist(lapply(seq_along(chromosomes), function(block) {
    chromosome <- chromosomes[[block]]
    vapply(seq_len(markers_per_chromosome), function(marker) {
      calls <- gt[((seq_along(donor_ids) + marker + block - 3L) %% 4L) + 1L]
      paste(c(chromosome, marker * 100L,
              paste0("chr", chromosome, "_m", marker),
              "A", "G", ".", "PASS", ".", "GT", calls), collapse = "\t")
    }, character(1L))
  }), use.names = FALSE)
  vcf <- file.path(root, "reference.vcf")
  writeLines(c("##fileformat=VCFv4.2", header, records), vcf, useBytes = TRUE)
  map <- data.frame(
    chromosome = rep(chromosomes, each = markers_per_chromosome),
    variant_id = variant_ids,
    genetic_position_cm = rep(seq(0, 150, length.out = markers_per_chromosome),
                              length(chromosomes)),
    stringsAsFactors = FALSE
  )

  reference <- gsim_import_vcf(
    vcf, map, file.path(root, "reference"), unsupported = "error"
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
  aligned_pedigree <- pedigree$pedigree[
    match(pedigree$canonical_order, pedigree$pedigree$animal), , drop = FALSE
  ]
  founder_ids <- aligned_pedigree$animal[
    is.na(aligned_pedigree$sire) & is.na(aligned_pedigree$dam)
  ]
  stopifnot(
    pedigree$diagnostics$full_sib_families > 0,
    pedigree$diagnostics$paternal_half_sib_sires > 0,
    pedigree$diagnostics$maternal_half_sib_dams > 0
  )

  populations <- setNames(rep(c("P1", "P2"), each = 4L), donor_ids)
  mutation_age <- setNames(rep(1e9, length(variant_ids)), variant_ids)
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

  Glist <- suppressMessages(qgg::gprep(
    study = "gsim-end-to-end",
    bedfiles = pedigree_bed$paths[["bed"]],
    bimfiles = pedigree_bed$paths[["bim"]],
    famfiles = pedigree_bed$paths[["fam"]]
  ))
  marker_ids <- as.character(unlist(Glist$rsids, use.names = FALSE))

  # gprep() supplies Glist$maf. This small positive named LD vector is only an
  # illustrative stand-in; realistic LD scores should be prepared once with
  # the genotype resource and aligned by marker ID.
  ld_score <- setNames(seq(0.8, 2.0, length.out = length(marker_ids)),
                       marker_ids)
  q <- setNames(seq(0.15, 0.75, length.out = length(marker_ids)), marker_ids)
  q[seq(1L, length(q), by = 8L)] <- 1
  annotation_score <- setNames(
    seq(0.75, 1.25, length.out = length(marker_ids)), marker_ids
  )

  default <- gsim(
    Glist = Glist, architecture = "bayesr", n_causal = 8L, seed = 8001
  )
  probability <- gsim(
    Glist = Glist, architecture = "bayesr",
    causal_probability = q, seed = 8002
  )
  variance <- gsim(
    Glist = Glist, architecture = "bayesr", n_causal = 8L,
    a = -0.4, b = -1, ld_score = ld_score, seed = 8001
  )
  combined <- gsim(
    Glist = Glist, architecture = "bayesr", causal_probability = q,
    a = -0.4, b = -1, c = 0.5,
    ld_score = ld_score, annotation_score = annotation_score,
    seed = 8002
  )

  stopifnot(
    identical(Glist$ids, pedigree$canonical_order),
    identical(rownames(default$Y), pedigree$canonical_order),
    identical(names(combined$causal_probability), marker_ids),
    identical(names(combined$marker_multipliers), marker_ids)
  )
  print(data.frame(
    model = c("default", "causal probability", "variance", "combined"),
    causal_markers = c(nrow(default$causal), nrow(probability$causal),
                       nrow(variance$causal), nrow(combined$causal))
  ))
  print(list(
    reference = reference$prefix,
    reusable_base = base$prefix,
    phased_pedigree = pedigree_hap_prefix,
    phenotype_Glist_samples = length(Glist$ids),
    phenotype_Glist_markers = length(marker_ids)
  ))
})
