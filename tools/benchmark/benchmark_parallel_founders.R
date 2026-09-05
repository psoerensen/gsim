#!/usr/bin/env Rscript

# Controlled founder-scaling diagnosis. The reference and all outputs must live
# outside the repository. One invocation measures one thread/batch setting.

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

env_integer <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, ".")
  }
  value
}
env_flag <- function(name, default = FALSE) {
  value <- tolower(Sys.getenv(name, if (default) "true" else "false"))
  if (!value %in% c("true", "false")) stop(name, " must be true or false.")
  identical(value, "true")
}

founder_count <- env_integer("GSIM_PARALLEL_FOUNDERS", 10000L)
batch_size <- env_integer("GSIM_PARALLEL_BATCH_SIZE", 4096L)
materialization_reps <- env_integer("GSIM_PARALLEL_MATERIALIZATION_REPS", 1L)
probe_reps <- env_integer("GSIM_PARALLEL_PROBE_REPS", 31L)
seed <- suppressWarnings(as.numeric(Sys.getenv(
  "GSIM_PARALLEL_SEED", "20260905"
)))
if (!is.finite(seed) || seed < 0 || batch_size %% 64L != 0L) {
  stop("Seed must be finite/nonnegative and batch size a multiple of 64.")
}
keep_output <- env_flag("GSIM_PARALLEL_KEEP_OUTPUT", FALSE)

timed <- function(expression) {
  before <- proc.time()
  value <- eval.parent(substitute(expression))
  delta <- proc.time() - before
  list(
    value = value,
    time = c(user = unname(delta[["user.self"]]),
             system = unname(delta[["sys.self"]]),
             elapsed = unname(delta[["elapsed"]]))
  )
}
add_time <- function(first, second) {
  stats::setNames(unname(first) + unname(second), names(first))
}
cpu_utilization <- function(value) {
  if (!is.finite(value[["elapsed"]]) || value[["elapsed"]] <= 0) return(NA_real_)
  100 * (value[["user"]] + value[["system"]]) / value[["elapsed"]]
}
coefficient_of_variation <- function(value) {
  if (length(value) < 2L || mean(value) == 0) return(0)
  stats::sd(value) / mean(value)
}
summarize_probe <- function(values) {
  c(minimum = min(values), median = stats::median(values), maximum = max(values))
}

plan_symbol <- get("C_gsim_hapnest_plan", envir = asNamespace("gsim"))
materialize_symbol <- get(
  "C_gsim_packed_materialize_founders", envir = asNamespace("gsim")
)

opened <- timed(gsim:::.gsim_hap_dataset_open(reference_prefix))
reader <- opened$value
on.exit(try(gsim:::.gsim_hap_dataset_close(reader), silent = TRUE), add = TRUE)
if (length(reader$chromosome) != 1L || reader$chromosome[[1L]] != "22") {
  stop("Benchmark reference must contain exactly chromosome 22.")
}
variants <- reader$variants
marker_count <- nrow(variants)
loaded <- timed(gsim:::.gsim_hap_dataset_load_chromosome(reader, "22"))
reference <- loaded$value
on.exit({
  try(gsim:::.gsim_packed_close(reference$h1), silent = TRUE)
  try(gsim:::.gsim_packed_close(reference$h2), silent = TRUE)
}, add = TRUE)

reference_ids <- reader$samples$individual_id
populations <- stats::setNames(rep.int("P1", length(reference_ids)), reference_ids)
mutation_age <- rep.int(1e9, marker_count)
founder_ids <- paste0("syn", seq_len(founder_count))
sample_metadata <- gsim:::.gsim_plink_sample_metadata(
  founder_ids, family_id = rep.int("base", founder_count)
)
output_prefix <- file.path(
  output_dir, paste0("founders-b", batch_size, "-t", threads)
)

sink <- gsim:::.gsim_hap_dataset_create(
  output_prefix, sample_metadata, overwrite = TRUE,
  provenance = list(operation = "founder scaling diagnosis",
                    threads = threads, batch_size = batch_size)
)
completed <- FALSE
on.exit(if (!completed) try(gsim:::.gsim_hap_dataset_cancel(sink), silent = TRUE),
        add = TRUE)
initialized <- timed(gsim:::.gsim_hap_dataset_begin_chromosome(
  sink, "22", variants
))

zero_time <- c(user = 0, system = 0, elapsed = 0)
plan_time <- zero_time
write_time <- zero_time
materialization_time <- matrix(
  0, nrow = materialization_reps, ncol = 3L,
  dimnames = list(paste0("rep", seq_len(materialization_reps)),
                  c("user", "system", "elapsed"))
)
max_event_bytes <- 0
max_batch_packed_bytes <- 0
starts <- seq.int(0L, founder_count - 1L, by = batch_size)
balance <- vector("list", length(starts))
first_plan <- NULL
batch_number <- 0L

workflow_start <- proc.time()
for (individual_offset in starts) {
  batch_number <- batch_number + 1L
  count <- min(batch_size, founder_count - individual_offset)
  input <- gsim:::.gsim_hapnest_packed_reference_inputs(
    reference$h1, reference$h2, populations, c(P1 = 1),
    c(P1 = length(reference_ids)), c(P1 = 10000), c(P1 = 0.02),
    variants$genetic_position_cm, mutation_age, count, seed, "22",
    return_genotypes = FALSE, return_segments = FALSE,
    individual_offset = individual_offset
  )
  planned <- timed(.Call(
    plan_symbol, input$donor_count, input$marker_count, input$donor_codes,
    input$weights, input$N, input$Ne, input$rho,
    rep.int(1L, input$marker_count), "22", input$genetic_position,
    input$n, input$seed, input$individual_offset
  ))
  plan <- planned$value
  plan_time <- add_time(plan_time, planned$time)
  max_event_bytes <- max(max_event_bytes, as.numeric(object.size(plan)))
  if (is.null(first_plan)) first_plan <- plan

  words <- ceiling(count / 64)
  worker_count <- min(threads, words)
  word_boundaries <- floor(words * (0:worker_count) / worker_count)
  local_individual <- as.numeric(plan$individual) - individual_offset - 1
  worker <- findInterval(
    local_individual, word_boundaries[-length(word_boundaries)] * 64
  )
  marker_visits <- as.numeric(plan$end) - as.numeric(plan$start) + 1
  balance[[batch_number]] <- data.frame(
    batch = batch_number,
    worker = seq_len(worker_count),
    output_words = diff(word_boundaries) * marker_count * 2,
    segment_operations = tabulate(worker, nbins = worker_count),
    marker_visits = as.numeric(rowsum(marker_visits, worker,
                                      reorder = FALSE)[, 1L]),
    stringsAsFactors = FALSE
  )

  h1 <- gsim:::.gsim_packed_zero(count, marker_count)
  h2 <- gsim:::.gsim_packed_zero(count, marker_count)
  for (repetition in seq_len(materialization_reps)) {
    materialized <- timed(.Call(
      materialize_symbol, h1, h2, reference$h1, reference$h2,
      plan$individual, plan$phase, plan$donor_individual, plan$start, plan$end,
      plan$coalescent_age, mutation_age, as.integer(individual_offset),
      as.integer(threads), FALSE
    ))
    materialization_time[repetition, ] <-
      materialization_time[repetition, ] + materialized$time
  }
  ids <- founder_ids[individual_offset + seq_len(count)]
  h1 <- gsim:::.gsim_packed_tag(h1, ids, variants$variant_id)
  h2 <- gsim:::.gsim_packed_tag(h2, ids, variants$variant_id)
  max_batch_packed_bytes <- max(
    max_batch_packed_bytes,
    unname(gsim:::.gsim_packed_info(h1)[[4L]] +
             gsim:::.gsim_packed_info(h2)[[4L]])
  )
  written <- timed(gsim:::.gsim_hap_dataset_write_batch(
    sink, h1, h2, individual_offset
  ))
  write_time <- add_time(write_time, written$time)
  gsim:::.gsim_packed_close(h1)
  gsim:::.gsim_packed_close(h2)
  rm(h1, h2, plan, input)
}

finalized <- timed(gsim:::.gsim_hap_dataset_finalize(sink))
manifest <- finalized$value
completed <- TRUE
workflow_delta <- proc.time() - workflow_start
workflow_measured <- c(
  user = unname(workflow_delta[["user.self"]]),
  system = unname(workflow_delta[["sys.self"]]),
  elapsed = unname(workflow_delta[["elapsed"]])
)
extra_repetitions <- if (materialization_reps > 1L) {
  colSums(materialization_time[-1L, , drop = FALSE])
} else zero_time
workflow_single_pass <- workflow_measured - extra_repetitions

# Bound serial validation and thread setup/join costs without changing native
# production code. The minimal plan touches one marker per worker.
probe_count <- min(batch_size, founder_count)
probe_h1 <- gsim:::.gsim_packed_zero(probe_count, marker_count)
probe_h2 <- gsim:::.gsim_packed_zero(probe_count, marker_count)
on.exit({
  try(gsim:::.gsim_packed_close(probe_h1), silent = TRUE)
  try(gsim:::.gsim_packed_close(probe_h2), silent = TRUE)
}, add = TRUE)
probe_words <- ceiling(probe_count / 64)
probe_workers <- min(threads, probe_words)
probe_destinations <- 1 + floor(
  probe_words * (0:(probe_workers - 1L)) / probe_workers
) * 64
probe_plan <- list(
  individual = as.numeric(probe_destinations),
  phase = rep.int(1L, probe_workers),
  donor_individual = rep.int(1L, probe_workers),
  start = rep.int(1L, probe_workers), end = rep.int(1L, probe_workers),
  coalescent_age = rep.int(0, probe_workers)
)
probe_call <- function(use_threads, plan) {
  timed(.Call(
    materialize_symbol, probe_h1, probe_h2, reference$h1, reference$h2,
    plan$individual, plan$phase, plan$donor_individual, plan$start, plan$end,
    plan$coalescent_age, mutation_age, 0L, as.integer(use_threads), FALSE
  ))$time[["elapsed"]]
}
probe_group <- function(use_threads, plan, repetitions) {
  timed(for (iteration in seq_len(repetitions)) {
    probe_call(use_threads, plan)
  })$time[["elapsed"]] / repetitions
}
empty_plan <- list(
  individual = numeric(), phase = integer(), donor_individual = integer(),
  start = integer(), end = integer(), coalescent_age = numeric()
)
validation_seconds <- replicate(5L, probe_group(1L, empty_plan, probe_reps))
minimal_one_seconds <- replicate(5L, probe_group(1L, probe_plan, probe_reps))
minimal_threaded_seconds <- replicate(
  5L, probe_group(threads, probe_plan, probe_reps)
)

# Preserve a realistic segment count while reducing every segment to one marker.
short_plan <- first_plan[c(
  "individual", "phase", "donor_individual", "coalescent_age"
)]
short_plan$start <- rep.int(1L, length(short_plan$phase))
short_plan$end <- short_plan$start
short_plan_seconds <- replicate(
  3L, probe_group(threads, short_plan, max(9L, probe_reps))
)

balance <- do.call(rbind, balance)
balance_totals <- stats::aggregate(
  cbind(output_words, segment_operations, marker_visits) ~ worker,
  data = balance, FUN = sum
)
reference_bytes <- unname(gsim:::.gsim_packed_info(reference$h1)[[4L]] +
                            gsim:::.gsim_packed_info(reference$h2)[[4L]])
paths <- manifest$paths
output_info <- list(
  paths = paths,
  bytes = unname(file.info(paths[["hap"]])$size),
  md5 = unname(tools::md5sum(paths)),
  individual_marker_combinations = founder_count * marker_count
)
median_materialization <- apply(materialization_time, 2L, stats::median)
result <- list(
  configuration = list(
    reference_prefix = reference_prefix, chromosome = "22",
    reference_individuals = length(reference_ids), markers = marker_count,
    founders = founder_count, batch_size = batch_size, threads = threads,
    batches = length(starts), materialization_repetitions = materialization_reps,
    probe_repetitions = probe_reps, N = length(reference_ids), Ne = 10000,
    rho = 0.02, ancestry_weights = c(P1 = 1), mutation_age = 1e9,
    seed = seed
  ),
  stage_times = list(
    hap_open = opened$time, hap_load = loaded$time,
    event_planning = plan_time,
    packed_materialization_repetitions = materialization_time,
    packed_materialization_median = median_materialization,
    hap_initialize = initialized$time, hap_batch_writes = write_time,
    hap_finalize = finalized$time,
    workflow_single_pass = workflow_single_pass,
    cpu_utilization_percent = c(
      materialization = cpu_utilization(median_materialization),
      workflow = cpu_utilization(workflow_single_pass)
    )
  ),
  worker_balance = list(
    per_batch = balance, totals = balance_totals,
    output_word_cv = coefficient_of_variation(balance_totals$output_words),
    segment_operation_cv = coefficient_of_variation(
      balance_totals$segment_operations
    ),
    marker_visit_cv = coefficient_of_variation(balance_totals$marker_visits),
    native_thread_creations_per_normal_workflow = sum(pmax(
      0, pmin(threads, ceiling(pmin(batch_size,
        founder_count - starts) / 64)) - 1L
    )),
    native_joins_per_normal_workflow = length(starts),
    partition_contract = paste(
      "contiguous disjoint 64-bit output words; reference planes shared",
      "read-only; no locks, atomics, R calls, or hot-loop allocation"
    )
  ),
  overhead_probes = list(
    zero_event_validation_seconds = summarize_probe(validation_seconds),
    minimal_one_thread_seconds = summarize_probe(minimal_one_seconds),
    minimal_requested_threads_seconds = summarize_probe(minimal_threaded_seconds),
    estimated_requested_thread_setup_join_seconds =
      stats::median(minimal_threaded_seconds) - stats::median(minimal_one_seconds),
    realistic_segment_count_one_marker_seconds = summarize_probe(
      short_plan_seconds
    )
  ),
  memory = list(
    reference_packed_bytes = reference_bytes,
    peak_batch_packed_bytes = max_batch_packed_bytes,
    maximum_event_plan_bytes = max_event_bytes,
    peak_reference_plus_batch_payload_bytes = reference_bytes +
      max_batch_packed_bytes,
    process_peak_working_set_bytes = suppressWarnings(as.numeric(Sys.getenv(
      "GSIM_BENCH_PEAK_WORKING_SET_BYTES", NA_character_
    )))
  ),
  output = output_info,
  rates = list(
    materialization_combinations_per_second = founder_count * marker_count /
      median_materialization[["elapsed"]],
    workflow_combinations_per_second = founder_count * marker_count /
      workflow_single_pass[["elapsed"]]
  ),
  session = utils::sessionInfo()
)

saveRDS(result, file.path(output_dir, "result.rds"))
dput(result[1:6], file = file.path(output_dir, "result.R"))
if (!keep_output) unlink(unname(paths), force = TRUE)
print(result[1:6])
