# Stage 6P: one fixed large pedigree and sequential record-view qualification.
# Run manually from the gsim repository after loading the development package.
# This script writes no output files and is never run by R CMD check.

if (!"package:gsim" %in% search()) library(gsim)

.format_bytes <- function(bytes) {
  format(structure(bytes, class = "object_size"), units = "auto")
}

.print_named <- function(label, value) {
  cat(label, ":\n", sep = "")
  print(value)
}

.pedigree_seed <- 600025L
.record_seeds <- c(single_trait = 600026L,
                   multitrait = 600027L,
                   longitudinal = 600028L)

cat("Stage 6P pedigree solver workload qualification\n")
cat("Pedigree seed:", .pedigree_seed, "\n")
cat("Record seeds:", paste(names(.record_seeds), .record_seeds,
                            sep = "=", collapse = ", "), "\n")

.pedigree_time <- system.time({
  pedigree <- gsim_pedigree(
    n_animals = 50000L,
    n_generations = 25L,
    founder_generations = 2L,
    sires_per_generation = 80L,
    dams_per_generation = 600L,
    overlapping_generation_probability = 0.20,
    unknown_sire_probability = 0.02,
    unknown_dam_probability = 0.02,
    new_founder_probability = 0.01,
    unphenotyped_probability = 0.15,
    permute_external_order = TRUE,
    seed = .pedigree_seed
  )
})

cat("\nPedigree summary\n")
.print_named("diagnostics", pedigree$diagnostics)
.print_named("checksums", pedigree$checksums)
cat("pedigree returned object size:",
    .format_bytes(as.numeric(object.size(pedigree))), "\n")
cat("pedigree construction elapsed seconds:", unname(.pedigree_time["elapsed"]),
    "\n")

.qualify_view <- function(model, seed) {
  elapsed <- system.time({
    view <- gsim_pedigree_records(
      pedigree,
      model = model,
      random_regression_rank = 3L,
      minimum_records = 1L,
      maximum_records = 5L,
      residual_groups = 3L,
      prediction_records = FALSE,
      permute_record_order = TRUE,
      seed = seed
    )
  })
  rank <- if (identical(model, "longitudinal")) 3L else
    if (identical(model, "multitrait")) 2L else 1L
  summary <- list(
    model = model,
    seed = seed,
    animals = pedigree$diagnostics$animals,
    founders = pedigree$diagnostics$founders,
    generations = pedigree$diagnostics$generations,
    sires_used = pedigree$diagnostics$sires_used,
    dams_used = pedigree$diagnostics$dams_used,
    unknown_sires = pedigree$diagnostics$unknown_sires,
    unknown_dams = pedigree$diagnostics$unknown_dams,
    unphenotyped = pedigree$diagnostics$unphenotyped,
    observation_units = view$diagnostics$observation_units,
    scalar_records = view$diagnostics$scalar_records,
    observation_patterns = view$diagnostics$observation_patterns,
    records_per_animal = summary(as.integer(view$diagnostics$records_per_animal)),
    residual_group_counts = view$diagnostics$residual_group_counts,
    fixed_design_dimensions = view$diagnostics$fixed_design$dimensions,
    fixed_design_stored_entries =
      view$diagnostics$fixed_design$stored_entries,
    basis_dimensions = view$diagnostics$basis_dimensions,
    basis_rank = view$diagnostics$basis_rank,
    approximate_mme_coefficient_dimension = view$fixed_design$ncol + rank *
      pedigree$diagnostics$animals,
    object_sizes_bytes = view$diagnostics$object_sizes_bytes,
    object_sizes_readable = vapply(view$diagnostics$object_sizes_bytes,
                                   .format_bytes, character(1L)),
    elapsed_seconds = unname(elapsed["elapsed"]),
    checksums = view$diagnostics$checksums,
    exact_relationship_draw = view$truth$exact_relationship_draw,
    purpose = view$truth$purpose
  )
  view <- NULL
  invisible(gc())
  summary
}

cat("\nSingle-trait view\n")
single_summary <- .qualify_view("single_trait", .record_seeds[["single_trait"]])
print(single_summary)

cat("\nTwo-trait incomplete-pattern view\n")
multitrait_summary <- .qualify_view("multitrait", .record_seeds[["multitrait"]])
print(multitrait_summary)

cat("\nLongitudinal rank-three view\n")
longitudinal_summary <- .qualify_view("longitudinal",
                                     .record_seeds[["longitudinal"]])
print(longitudinal_summary)

cat("\nQualification complete. No MME, relationship matrix, REML/BLUP fit,",
    "factorization, or solver benchmark was constructed.\n")
