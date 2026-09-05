#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("usage: benchmark_parallel_founders.R REFERENCE_PREFIX OUTPUT_DIR THREADS")
}
reference_prefix <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
output_dir <- normalizePath(args[[2L]], winslash = "/", mustWork = FALSE)
threads <- suppressWarnings(as.integer(args[[3L]]))
if (is.na(threads) || threads < 1L) stop("THREADS must be a positive integer.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
repository <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (startsWith(paste0(output_dir, "/"), paste0(repository, "/"))) {
  stop("Benchmark output must be outside the repository.")
}
suppressPackageStartupMessages(library(gsim))
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1")

founder_count <- as.integer(Sys.getenv("GSIM_PARALLEL_FOUNDERS", "10000"))
batch_size <- as.integer(Sys.getenv("GSIM_PARALLEL_BATCH_SIZE", "4096"))
seed <- as.numeric(Sys.getenv("GSIM_PARALLEL_SEED", "20260905"))
if (is.na(founder_count) || founder_count < 1L || is.na(batch_size) ||
    batch_size < 1L || batch_size %% 64L != 0L) {
  stop("Founder count must be positive and batch size a positive multiple of 64.")
}

elapsed <- function(expression) {
  before <- proc.time()
  value <- eval.parent(substitute(expression))
  list(value = value, time = proc.time() - before)
}
seconds <- function(x) unname(x[["elapsed"]])

backend <- gsim:::.gsim_packed_backend()
metadata_backend <- gsim:::.gsim_metadata_backend()
total_start <- proc.time()
opened <- elapsed(gsim:::.gsim_hap_dataset_open(
  backend, metadata_backend, reference_prefix
))
reader <- opened$value
on.exit(try(gsim:::.gsim_hap_dataset_close(reader), silent = TRUE), add = TRUE)
if (length(reader$chromosome) != 1L || reader$chromosome[[1L]] != "22") {
  stop("Benchmark reference must contain exactly chromosome 22.")
}
variants <- reader$variants
marker_count <- nrow(variants)
loaded <- elapsed(gsim:::.gsim_hap_dataset_load_chromosome(reader, "22"))
reference <- loaded$value
on.exit({
  try(gsim:::.gsim_packed_close(reference$h1), silent = TRUE)
  try(gsim:::.gsim_packed_close(reference$h2), silent = TRUE)
}, add = TRUE)
reference_ids <- reader$samples$individual_id
populations <- rep.int("P1", length(reference_ids))
names(populations) <- reference_ids
mutation_age <- rep.int(1e9, marker_count)
founder_ids <- paste0("syn", seq_len(founder_count))
sample_metadata <- gsim:::.gsim_plink_sample_metadata(
  founder_ids, family_id = rep.int("base", founder_count)
)
output_prefix <- file.path(output_dir, paste0("founders-t", threads))
sink <- gsim:::.gsim_hap_dataset_create(
  backend, metadata_backend, output_prefix, sample_metadata,
  overwrite = TRUE,
  provenance = list(operation = "parallel founder benchmark", threads = threads)
)
completed <- FALSE
on.exit(if (!completed) try(gsim:::.gsim_hap_dataset_cancel(sink), silent = TRUE),
        add = TRUE)
initialized <- elapsed(gsim:::.gsim_hap_dataset_begin_chromosome(
  sink, "22", variants
))
plan_seconds <- materialization_seconds <- write_seconds <- 0
max_event_bytes <- 0
max_batch_packed_bytes <- 0
plan_symbol <- get("C_gsim_hapnest_plan", envir = asNamespace("gsim"))
materialize_symbol <- get("C_gsim_packed_materialize_founders",
                          envir = asNamespace("gsim"))
starts <- seq.int(0L, founder_count - 1L, by = batch_size)
for (individual_offset in starts) {
  count <- min(batch_size, founder_count - individual_offset)
  input <- gsim:::.gsim_hapnest_packed_reference_inputs(
    backend, reference$h1, reference$h2, populations, c(P1 = 1),
    c(P1 = length(reference_ids)), c(P1 = 10000), c(P1 = 0.02),
    variants$genetic_position_cm, mutation_age, count, seed, "22",
    return_genotypes = FALSE, return_segments = FALSE,
    individual_offset = individual_offset
  )
  planned <- elapsed(.Call(
    plan_symbol, input$donor_count, input$marker_count, input$donor_codes,
    input$weights, input$N, input$Ne, input$rho,
    rep.int(1L, input$marker_count), "22", input$genetic_position,
    input$n, input$seed, input$individual_offset
  ))
  plan <- planned$value
  plan_seconds <- plan_seconds + seconds(planned$time)
  max_event_bytes <- max(max_event_bytes, as.numeric(object.size(plan)))
  h1 <- gsim:::.gsim_packed_zero(backend, count, marker_count)
  h2 <- gsim:::.gsim_packed_zero(backend, count, marker_count)
  materialized <- elapsed(.Call(
    materialize_symbol, h1, h2, reference$h1, reference$h2,
    plan$individual, plan$phase, plan$donor_individual, plan$start, plan$end,
    plan$coalescent_age, mutation_age, as.integer(individual_offset),
    as.integer(threads), FALSE
  ))
  materialization_seconds <- materialization_seconds + seconds(materialized$time)
  ids <- founder_ids[individual_offset + seq_len(count)]
  h1 <- gsim:::.gsim_packed_tag(h1, ids, variants$variant_id)
  h2 <- gsim:::.gsim_packed_tag(h2, ids, variants$variant_id)
  max_batch_packed_bytes <- max(
    max_batch_packed_bytes,
    unname(gsim:::.gsim_packed_info(h1)[[4L]] +
             gsim:::.gsim_packed_info(h2)[[4L]])
  )
  written <- elapsed(gsim:::.gsim_hap_dataset_write_batch(
    sink, h1, h2, individual_offset
  ))
  write_seconds <- write_seconds + seconds(written$time)
  gsim:::.gsim_packed_close(h1)
  gsim:::.gsim_packed_close(h2)
  rm(h1, h2, plan, input)
}
finalized <- elapsed(gsim:::.gsim_hap_dataset_finalize(sink))
manifest <- finalized$value
completed <- TRUE
total <- proc.time() - total_start
reference_bytes <- unname(gsim:::.gsim_packed_info(reference$h1)[[4L]] +
                            gsim:::.gsim_packed_info(reference$h2)[[4L]])
result <- list(
  configuration = list(
    reference_prefix = reference_prefix, chromosome = "22",
    reference_individuals = length(reference_ids), markers = marker_count,
    founders = founder_count, batch_size = batch_size, threads = threads,
    N = length(reference_ids), Ne = 10000, rho = 0.02,
    ancestry_weights = c(P1 = 1), mutation_age = 1e9, seed = seed
  ),
  timings = c(
    hap_open = seconds(opened$time), hap_load = seconds(loaded$time),
    event_planning = plan_seconds,
    packed_materialization = materialization_seconds,
    hap_initialize = seconds(initialized$time), hap_batch_writes = write_seconds,
    hap_finalize_metadata_publish = seconds(finalized$time),
    total = seconds(total)
  ),
  memory = list(
    reference_packed_bytes = reference_bytes,
    peak_batch_packed_bytes = max_batch_packed_bytes,
    maximum_event_plan_bytes = max_event_bytes,
    peak_reference_plus_batch_payload_bytes = reference_bytes +
      max_batch_packed_bytes,
    process_peak_working_set_bytes = as.numeric(Sys.getenv(
      "GSIM_BENCH_PEAK_WORKING_SET_BYTES", NA_character_))
  ),
  output = list(
    path = manifest$paths[["hap"]], bytes = file.info(manifest$paths[["hap"]])$size,
    md5 = unname(tools::md5sum(manifest$paths[["hap"]])),
    individual_marker_updates = founder_count * marker_count,
    updates_per_second = founder_count * marker_count / seconds(total)
  ),
  session = utils::sessionInfo()
)
saveRDS(result, file.path(output_dir, paste0("result-t", threads, ".rds")))
dput(result[1:4], file = file.path(output_dir,
                                   paste0("result-t", threads, ".R")))
print(result[1:4])
