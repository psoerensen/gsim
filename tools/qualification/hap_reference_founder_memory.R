# Bounded, reproducible accounting for the direct HAP-reference founder path.
# Requires GSIM_GBITS_LIBRARY and GSIM_GMAT_LIBRARY; this is not a benchmark.

local({
library(gsim)

backend <- gsim:::.gsim_gbits_backend(Sys.getenv("GSIM_GBITS_LIBRARY"))
metadata_backend <- gsim:::.gsim_gmat_backend(Sys.getenv("GSIM_GMAT_LIBRARY"))
root <- tempfile("gsim-hap-reference-memory-", tmpdir = Sys.getenv("TEMP"))
dir.create(root)
on.exit(unlink(root, recursive = TRUE), add = TRUE)

donors <- paste0("r", seq_len(128L))
chromosomes <- c(`Z-2` = 129L, `01` = 65L, chrA = 33L)
prefix <- file.path(root, "reference")
dataset <- gsim:::.gsim_hap_dataset_create(
  backend, metadata_backend, prefix,
  gsim:::.gsim_plink_sample_metadata(donors))
for (block in seq_along(chromosomes)) {
  label <- names(chromosomes)[[block]]
  marker_count <- chromosomes[[block]]
  variant_ids <- paste0(label, "_v", seq_len(marker_count))
  h1 <- outer(seq_along(donors), seq_len(marker_count), function(i, j) {
    (i * 3L + j * 5L + block) %% 2L
  })
  h2 <- outer(seq_along(donors), seq_len(marker_count), function(i, j) {
    (i * 7L + j * 3L + j %/% 2L + block + 1L) %% 2L
  })
  dimnames(h1) <- dimnames(h2) <- list(donors, variant_ids)
  packed_h1 <- gsim:::.gsim_gbits_pack(backend, h1)
  packed_h2 <- gsim:::.gsim_gbits_pack(backend, h2)
  position <- seq(0, 2, length.out = marker_count)
  variants <- data.frame(
    chromosome = rep.int(label, marker_count), variant_id = variant_ids,
    genetic_position_cm = position,
    base_pair_position = seq_len(marker_count),
    alt = rep.int("C", marker_count), ref = rep.int("A", marker_count),
    bit1_allele = rep.int("C", marker_count),
    bit0_allele = rep.int("A", marker_count), stringsAsFactors = FALSE)
  gsim:::.gsim_hap_dataset_append(
    dataset, label, packed_h1, packed_h2, variants)
  gsim:::.gsim_gbits_close(packed_h1)
  gsim:::.gsim_gbits_close(packed_h2)
}
gsim:::.gsim_hap_dataset_finalize(dataset)
rm(h1, h2)
gc()

reader <- gsim:::.gsim_hap_dataset_open(backend, metadata_backend, prefix)
population <- stats::setNames(rep(c("A", "B"), each = 64L), donors)
rows <- vector("list", length(chromosomes))
for (block in seq_along(chromosomes)) {
  label <- names(chromosomes)[[block]]
  marker_count <- chromosomes[[block]]
  variant_ids <- paste0(label, "_v", seq_len(marker_count))
  range <- reader$variants$chromosome == label
  position <- stats::setNames(
    reader$variants$genetic_position_cm[range], variant_ids)
  mutation_age <- stats::setNames(
    rep(c(0.25, 2, 20, 1e9), length.out = marker_count), variant_ids)
  out <- gsim:::.gsim_hapnest_founders_from_hap_chromosome(
    reader, label, population, c(A = 0.4, B = 0.6),
    c(A = 64, B = 64), c(A = 400, B = 700), c(A = 0.8, B = 1.3),
    position, mutation_age, n = 128L, seed = 20260905L,
    return_genotypes = FALSE, return_segments = TRUE)
  rows[[block]] <- data.frame(
    chromosome = label, markers = marker_count,
    reference_packed_bytes = out$memory$reference_packed_bytes,
    generated_packed_bytes = out$memory$generated_packed_bytes,
    event_record_bytes = out$memory$event_record_bytes,
    peak_biological_payload_bytes = out$memory$peak_biological_payload_bytes,
    raw_reference_bytes_avoided = out$memory$reference_raw_bytes_avoided,
    dense_genotype_bytes_avoided = 128 * marker_count,
    segment_count = nrow(out$segments), stringsAsFactors = FALSE)
  gsim:::.gsim_gbits_close(out$h1)
  gsim:::.gsim_gbits_close(out$h2)
}
gsim:::.gsim_hap_dataset_close(reader)
print(do.call(rbind, rows), row.names = FALSE)
})
