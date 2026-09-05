#!/usr/bin/env Rscript

# Reproducible chromosome 22 production benchmark for standalone gsim.
# Large inputs and outputs must live outside the repository. Downloads occur
# only when GSIM_BENCH_DOWNLOAD=true is set explicitly.

options(stringsAsFactors = FALSE)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1")

truthy <- function(x) tolower(x) %in% c("1", "true", "yes")
env_integer <- function(name, default) {
  value <- Sys.getenv(name, as.character(default))
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || parsed < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  parsed
}
env_number <- function(name, default) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || !is.finite(value)) {
    stop(name, " must be finite numeric.", call. = FALSE)
  }
  value
}

vcf_url <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
  "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.",
  "20130502.genotypes.vcf.gz"
)
map_url <- paste0(
  "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/",
  "20110106_recombination_hotspots/",
  "HapmapII_GRCh37_RecombinationHotspots.tar.gz"
)
data_dir <- Sys.getenv("GSIM_BENCH_DATA_DIR")
output_dir <- Sys.getenv("GSIM_BENCH_OUTPUT_DIR")
if (!nzchar(data_dir) || !nzchar(output_dir)) {
  stop("Set GSIM_BENCH_DATA_DIR and GSIM_BENCH_OUTPUT_DIR outside the repository.",
       call. = FALSE)
}
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
repository <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (startsWith(paste0(output_dir, "/"), paste0(repository, "/"))) {
  stop("Benchmark output directory must be outside the repository.", call. = FALSE)
}

vcf <- file.path(data_dir, basename(vcf_url))
map_archive <- file.path(data_dir, basename(map_url))
map_file <- file.path(data_dir, "genetic_map_GRCh37_chr22.txt")
allow_download <- truthy(Sys.getenv("GSIM_BENCH_DOWNLOAD", "false"))
if (!file.exists(vcf) && allow_download) {
  download.file(vcf_url, vcf, mode = "wb", quiet = FALSE)
}
if (!file.exists(map_archive) && !file.exists(map_file) && allow_download) {
  download.file(map_url, map_archive, mode = "wb", quiet = FALSE)
}
if (!file.exists(map_file) && file.exists(map_archive)) {
  utils::untar(map_archive, files = basename(map_file), exdir = data_dir)
}
if (!file.exists(vcf) || !file.exists(map_file)) {
  stop("Official VCF and extracted GRCh37 chromosome 22 map are required; ",
       "set GSIM_BENCH_DOWNLOAD=true for the explicitly authorized download.",
       call. = FALSE)
}

suppressPackageStartupMessages(library(gsim))
reference_count <- env_integer("GSIM_BENCH_REFERENCE_COUNT", 500L)
founder_count <- env_integer("GSIM_BENCH_FOUNDER_COUNT", 1000L)
region_start <- env_integer("GSIM_BENCH_REGION_START", 20000000L)
region_end <- env_integer("GSIM_BENCH_REGION_END", 22200000L)
if (region_start > region_end) stop("Benchmark region is invalid.", call. = FALSE)
seed <- env_number("GSIM_BENCH_SEED", 20260905)
rho <- env_number("GSIM_BENCH_RHO", 0.02)
Ne <- env_number("GSIM_BENCH_NE", 10000)
mutation_value <- env_number("GSIM_BENCH_MUTATION_AGE", 1e9)
calibration <- truthy(Sys.getenv("GSIM_BENCH_CALIBRATION", "false"))

# Header-only setup is not measured. It determines an exact, reproducible panel.
connection <- gzfile(vcf, "rt")
on.exit(try(close(connection), silent = TRUE), add = TRUE)
repeat {
  line <- readLines(connection, n = 1L, warn = FALSE)
  if (!length(line)) stop("VCF #CHROM header was not found.", call. = FALSE)
  if (startsWith(line, "#CHROM\t")) break
}
close(connection)
header <- strsplit(line, "\t", fixed = TRUE)[[1L]]
all_samples <- header[-seq_len(9L)]
if (reference_count > length(all_samples)) {
  stop("Requested more reference samples than the VCF contains.", call. = FALSE)
}
selected_samples <- all_samples[seq_len(reference_count)]
writeLines(selected_samples, file.path(output_dir, "selected_samples.txt"),
           useBytes = TRUE)

map_source <- utils::read.delim(map_file, check.names = FALSE)
if (!identical(unique(map_source[["Chromosome"]]), "chr22")) {
  stop("The authoritative map does not have the expected chr22 label.", call. = FALSE)
}
# Deliberate caller-side conversion; gsim itself performs exact label matching.
genetic_map <- data.frame(
  chromosome = "22",
  base_pair_position = map_source[["Position(bp)"]],
  genetic_position_cm = map_source[["Map(cM)"]]
)

timings <- list()
current_rss <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  unname(ps::ps_memory_info(ps::ps_handle())[["rss"]])
}
stage <- function(name, expression) {
  gc()
  before <- proc.time()
  value <- eval.parent(substitute(expression))
  elapsed <- proc.time() - before
  timings[[name]] <<- c(
    user_seconds = unname(elapsed[["user.self"]]),
    system_seconds = unname(elapsed[["sys.self"]]),
    wall_seconds = unname(elapsed[["elapsed"]]),
    rss_after_bytes = current_rss()
  )
  value
}

# Warm-up: native symbol loading and a tiny packed round trip only. It does not
# scan the real VCF and is excluded from all reported timings.
packed_backend <- gsim:::.gsim_packed_backend()
metadata_backend <- gsim:::.gsim_metadata_backend()
warm <- matrix(as.raw(c(0, 1, 1, 0)), 2L, 2L)
warm_handle <- gsim:::.gsim_packed_pack(packed_backend, warm)
stopifnot(identical(gsim:::.gsim_packed_unpack(warm_handle), warm))
gsim:::.gsim_packed_close(warm_handle)
gc()

reference_prefix <- file.path(output_dir, "reference")
reference <- stage("vcf_to_hap_total", gsim_import_vcf(
  vcf = vcf, map = genetic_map, output = reference_prefix,
  samples = selected_samples, chromosome = "22",
  region = c(region_start, region_end), unsupported = "skip",
  overwrite = TRUE
))
reader <- stage("prepared_hap_open", gsim:::.gsim_hap_dataset_open(
  packed_backend, metadata_backend, reference_prefix
))
on.exit(try(gsim:::.gsim_hap_dataset_close(reader), silent = TRUE), add = TRUE)
if (!identical(reader$samples$individual_id, selected_samples)) {
  stop("Prepared HAP sample order differs from the frozen panel.", call. = FALSE)
}
if (!identical(reader$chromosome, "22")) {
  stop("Prepared HAP chromosome identity differs from the request.", call. = FALSE)
}
marker_count <- nrow(reader$variants)
if (!calibration && (marker_count < 25000L || marker_count > 100000L)) {
  stop("Frozen region retained ", marker_count,
       " markers, outside the predeclared 25,000-100,000 bound.", call. = FALSE)
}
variant_ids <- reader$variants$variant_id
positions_cm <- stats::setNames(reader$variants$genetic_position_cm, variant_ids)
mutation_age <- stats::setNames(rep(mutation_value, marker_count), variant_ids)
populations <- stats::setNames(rep("P1", reference_count), selected_samples)

reference_handles <- stage("prepared_hap_load_chromosome",
  gsim:::.gsim_hap_dataset_load_chromosome(reader, "22"))
on.exit({
  try(gsim:::.gsim_packed_close(reference_handles$h1), silent = TRUE)
  try(gsim:::.gsim_packed_close(reference_handles$h2), silent = TRUE)
}, add = TRUE)

founder_input <- gsim:::.gsim_hapnest_packed_reference_inputs(
  packed_backend, reference_handles$h1, reference_handles$h2,
  unname(populations), c(P1 = 1), c(P1 = reference_count), c(P1 = Ne),
  c(P1 = rho), unname(positions_cm), unname(mutation_age), founder_count,
  seed, rep.int("22", marker_count), "hapnest", FALSE, TRUE, 0L
)
plan_symbol <- get("C_gsim_hapnest_plan", envir = asNamespace("gsim"))
event_plan <- stage("founder_event_plan_only", .Call(
  plan_symbol, founder_input$donor_count, founder_input$marker_count,
  founder_input$donor_codes, founder_input$weights, founder_input$N,
  founder_input$Ne, founder_input$rho, rep.int(1L, marker_count), "22",
  founder_input$genetic_position, founder_input$n, founder_input$seed,
  founder_input$individual_offset
))

founders <- stage("founder_plan_and_packed_materialization",
  gsim:::.gsim_hapnest_founders_packed_reference_chromosome(
    packed_backend, reference_handles$h1, reference_handles$h2,
    unname(populations), c(P1 = 1), c(P1 = reference_count), c(P1 = Ne),
    c(P1 = rho), unname(positions_cm), unname(mutation_age), founder_count,
    seed, "22", return_genotypes = FALSE, return_segments = TRUE
  ))
on.exit({
  try(gsim:::.gsim_packed_close(founders$h1), silent = TRUE)
  try(gsim:::.gsim_packed_close(founders$h2), silent = TRUE)
}, add = TRUE)
if (!identical(event_plan$individual, founders$segments$individual) ||
    !identical(event_plan$phase, founders$segments$phase) ||
    !identical(event_plan$start, founders$segments$start) ||
    !identical(event_plan$end, founders$segments$end) ||
    !identical(event_plan$donor_individual, founders$segments$donor_individual) ||
    !identical(event_plan$coalescent_age, founders$segments$coalescent_age) ||
    !identical(event_plan$sampled_length, founders$segments$sampled_length)) {
  stop("Standalone event plan and materialized founder audit differ.", call. = FALSE)
}

# Reference handles can be released as soon as founders are materialized.
reference_info <- gsim:::.gsim_packed_info(reference_handles$h1)
gsim:::.gsim_packed_close(reference_handles$h1)
gsim:::.gsim_packed_close(reference_handles$h2)
reference_handles <- list(h1 = NULL, h2 = NULL)

founder_ids <- founders$sample_ids
sample_metadata <- gsim:::.gsim_plink_sample_metadata(
  founder_ids, family_id = rep.int("founders", founder_count)
)
variants <- reader$variants
founder_hap_prefix <- file.path(output_dir, "founders-hap")
founder_hap <- stage("founder_hap_output", {
  sink <- gsim:::.gsim_hap_dataset_create(
    packed_backend, metadata_backend, founder_hap_prefix, sample_metadata,
    overwrite = TRUE, provenance = list(benchmark = "1000G chr22"))
  gsim:::.gsim_hap_dataset_append(sink, "22", founders$h1, founders$h2,
                                  variants)
  gsim:::.gsim_hap_dataset_finalize(sink)
})
founder_bed_prefix <- file.path(output_dir, "founders-bed")
founder_bed <- stage("founder_bed_output", {
  sink <- gsim:::.gsim_plink_dataset_create(
    packed_backend, metadata_backend, founder_bed_prefix, sample_metadata,
    overwrite = TRUE, provenance = list(benchmark = "1000G chr22"))
  gsim:::.gsim_plink_dataset_append(sink, "22", founders$h1, founders$h2,
                                    variants)
  gsim:::.gsim_plink_dataset_finalize(sink)
})

pedigree <- stage("pedigree_structure", gsim_pedigree(
  n_generations = 2L, animals_per_generation = founder_count,
  founder_generations = 1L, unknown_sire_probability = 0,
  unknown_dam_probability = 0, new_founder_probability = 0,
  permute_external_order = TRUE, seed = 771
))
pedigree_founders <- pedigree$canonical_order[seq_len(founder_count)]
founders$h1 <- gsim:::.gsim_packed_tag(founders$h1, pedigree_founders, variant_ids)
founders$h2 <- gsim:::.gsim_packed_tag(founders$h2, pedigree_founders, variant_ids)
descendants <- stage("packed_pedigree_meiosis",
  gsim:::.gsim_pedigree_genotypes_packed_chromosome(
    packed_backend, pedigree, list(h1 = founders$h1, h2 = founders$h2),
    rep.int("22", marker_count), variants$genetic_position_cm / 100, seed,
    return_haplotypes = TRUE, return_genotypes = FALSE,
    return_crossovers = TRUE
  ))
on.exit({
  try(gsim:::.gsim_packed_close(descendants$h1), silent = TRUE)
  try(gsim:::.gsim_packed_close(descendants$h2), silent = TRUE)
}, add = TRUE)

repeat_reference <- stage("repeat_hap_load_chromosome",
  gsim:::.gsim_hap_dataset_load_chromosome(reader, "22"))
repeat_founders <- stage("repeat_founder_simulation",
  gsim:::.gsim_hapnest_founders_packed_reference_chromosome(
    packed_backend, repeat_reference$h1, repeat_reference$h2,
    unname(populations), c(P1 = 1), c(P1 = reference_count), c(P1 = Ne),
    c(P1 = rho), unname(positions_cm), unname(mutation_age), founder_count,
    seed + 1, "22", return_genotypes = FALSE, return_segments = TRUE
  ))
gsim:::.gsim_packed_close(repeat_reference$h1)
gsim:::.gsim_packed_close(repeat_reference$h2)
gsim:::.gsim_packed_close(repeat_founders$h1)
gsim:::.gsim_packed_close(repeat_founders$h2)
repeat_reference <- repeat_founders <- NULL

packed_summary <- function(h1, h2, individuals, markers) {
  chosen <- unique(as.integer(round(seq(1, marker_count,
                                         length.out = min(markers, marker_count)))))
  words <- ceiling(individuals / 64)
  phase <- c(h1 = 0, h2 = 0)
  dosage <- integer(3L)
  observed <- 0L
  for (marker in chosen) {
    for (word in seq_len(words)) {
      first <- (word - 1L) * 64L + 1L
      valid <- min(64L, individuals - first + 1L)
      b1 <- as.integer(rawToBits(gsim:::.gsim_packed_word(h1, marker, word)))[seq_len(valid)]
      b2 <- as.integer(rawToBits(gsim:::.gsim_packed_word(h2, marker, word)))[seq_len(valid)]
      phase <- phase + c(sum(b1), sum(b2))
      dosage <- dosage + tabulate(b1 + b2 + 1L, nbins = 3L)
      observed <- observed + valid
    }
  }
  list(markers = chosen, phase_allele_frequency = phase / observed,
       dosage_frequency = stats::setNames(dosage / observed, 0:2))
}
founder_summary <- packed_summary(founders$h1, founders$h2,
                                  founder_count, 256L)
segment_count <- table(interaction(founders$segments$individual,
                                   founders$segments$phase, drop = TRUE))
copied <- sum(founders$segments$copied_alternative)
retained <- sum(founders$segments$retained_alternative)

timing_table <- do.call(rbind, timings)
configuration <- list(
  region = c(region_start, region_end), chromosome = "22",
  reference_individuals = reference_count, reference_samples = selected_samples,
  reference_sample_file_md5 = unname(tools::md5sum(
    file.path(output_dir, "selected_samples.txt"))),
  marker_count = marker_count, variant_first = variant_ids[[1L]],
  variant_last = variant_ids[[marker_count]], founders = founder_count,
  pedigree_individuals = length(pedigree$canonical_order),
  population = "P1", ancestry_weights = c(P1 = 1), N = reference_count,
  Ne = Ne, rho = rho, mutation_age = mutation_value, seed = seed,
  threads = 1L, map_build = "GRCh37",
  map_label_conversion = "caller-side chr22 to exact VCF label 22"
)
configuration$calibration <- calibration
memory <- list(
  reference_packed_bytes = 2 * 8 * marker_count * ceiling(reference_count / 64),
  founder_packed_bytes = 2 * 8 * marker_count * ceiling(founder_count / 64),
  pedigree_packed_bytes = 2 * 8 * marker_count *
    ceiling(length(pedigree$canonical_order) / 64),
  event_plan_bytes = as.numeric(object.size(event_plan)),
  founder_audit_bytes = as.numeric(object.size(founders$segments)),
  import_peak_chromosome_packed_bytes =
    reference$import$peak_chromosome_packed_bytes,
  process_peak_working_set_bytes = suppressWarnings(as.numeric(
    Sys.getenv("GSIM_BENCH_PEAK_WORKING_SET_BYTES", NA_character_)))
)
files <- list(
  input_vcf_bytes = unname(file.info(vcf)$size),
  input_vcf_md5 = unname(tools::md5sum(vcf)),
  map_file_bytes = unname(file.info(map_file)$size),
  map_file_md5 = unname(tools::md5sum(map_file)),
  reference = stats::setNames(unname(file.info(reference$paths)$size),
                              names(reference$paths)),
  founder_hap = stats::setNames(unname(file.info(founder_hap$paths)$size),
                                names(founder_hap$paths)),
  founder_bed = stats::setNames(unname(file.info(founder_bed$paths)$size),
                                names(founder_bed$paths))
)
sanity <- list(
  two_phases = TRUE,
  phase_specific_allele_frequency = founder_summary$phase_allele_frequency,
  genotype_dosage_frequency = founder_summary$dosage_frequency,
  sampled_marker_count = length(founder_summary$markers),
  simulated_haplotypes = length(segment_count),
  segment_records = nrow(founders$segments),
  segment_count_per_haplotype = summary(as.numeric(segment_count)),
  sampled_length = summary(founders$segments$sampled_length),
  copied_genetic_span = summary(founders$segments$copied_genetic_span),
  donor_population_proportion = prop.table(table(founders$segments$donor_population)),
  mutation_retention_fraction = if (copied == 0) NA_real_ else retained / copied,
  pedigree_crossover_records = nrow(descendants$crossover_audit$crossovers),
  dense_genotype_matrix_allocated = FALSE
)
hapnest <- list(
  requested_revision = "ba52da1a63cf609306ea92540b3d130fa1efd213",
  source_directory = Sys.getenv("GSIM_BENCH_HAPNEST_DIR", ""),
  julia = Sys.which("julia"),
  executed = FALSE,
  blocker = if (!nzchar(Sys.which("julia"))) "Julia executable unavailable" else
    "Pinned local HAPNEST checkout/prepared inputs unavailable"
)
result <- list(
  configuration = configuration, timings = timing_table,
  import_report = reference$import, memory = memory, files = files,
  sanity = sanity, hapnest = hapnest,
  session = list(R = R.version.string, platform = R.version$platform,
                 os = Sys.info(), packages = utils::sessionInfo(),
                 environment_threads = Sys.getenv(
                   c("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")))
)
saveRDS(result, file.path(output_dir, "benchmark_result.rds"))
dput(result, file = file.path(output_dir, "benchmark_result.R"))
write.csv(timing_table, file.path(output_dir, "timings.csv"), quote = FALSE)
print(configuration[c("region", "reference_individuals", "marker_count",
                      "founders", "pedigree_individuals", "threads")])
print(timing_table)
print(memory)
print(files)
print(sanity)
