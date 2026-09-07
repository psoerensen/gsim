# BED-backed Glist integration. Scientific transformations stay in gsim's R
# routines; only the selected packed scalar product is native.
.gsim_bed_plan <- function(Glist, marker_ids, sample_ids) {
  unique_ids <- function(x, label) {
    x <- as.character(x)
    if (!length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x))
      .gsim_stop(label, " must be nonempty, unique, nonmissing IDs.")
    x
  }
  available <- unique_ids(Glist$ids, "Glist$ids")
  sample_ids <- unique_ids(sample_ids, "Selected sample IDs")
  if (anyNA(match(sample_ids, available))) .gsim_stop("Selected samples are absent from Glist$ids.")
  if (!is.null(Glist$n) && (!is.numeric(Glist$n) || length(Glist$n) != 1L ||
      is.na(Glist$n) || Glist$n != length(available))) .gsim_stop("Glist$n must match Glist$ids.")
  beds <- Glist$bedfiles
  if (!is.character(beds) || !length(beds) || anyNA(beds) || any(!nzchar(beds)))
    .gsim_stop("Native Glist simulation requires SNP-major BED files in Glist$bedfiles; ",
               "other storage requires an explicit deterministic getG_fun.")
  if (!is.list(Glist$rsids) || length(Glist$rsids) != length(beds))
    .gsim_stop("Glist$rsids must be a list aligned with Glist$bedfiles.")
  catalog <- unique_ids(unlist(Glist$rsids, use.names = FALSE), "Glist$rsids")
  if (anyNA(match(marker_ids, catalog))) .gsim_stop("Simulation markers are absent from Glist$rsids.")
  companion <- function(field, suffix) {
    value <- Glist[[field]]
    if (is.null(value)) {
      if (any(!grepl("[.]bed$", beds))) .gsim_stop("Supply Glist$", field, " for BED paths without .bed suffix.")
      value <- sub("[.]bed$", suffix, beds)
    }
    if (!is.character(value) || length(value) != length(beds) || anyNA(value) || any(!nzchar(value)))
      .gsim_stop("Glist$", field, " must align with bedfiles.")
    value
  }
  bims <- companion("bimfiles", ".bim"); fams <- companion("famfiles", ".fam")
  read_metadata <- function(path) {
    # Simulation FAM validators enforce parent-before-child pedigrees, which
    # ordinary qgg resources need not have. Only physical identity is needed here.
    tab <- utils::read.table(path, header = FALSE, colClasses = "character",
                             quote = "", comment.char = "", stringsAsFactors = FALSE)
    if (ncol(tab) != 6L || !nrow(tab)) .gsim_stop("Expected six-column PLINK metadata: ", path)
    tab
  }
  files <- vector("list", length(beds))
  file_index <- variant_index <- integer(length(marker_ids))
  for (f in seq_along(beds)) {
    fam <- read_metadata(fams[f]); physical_ids <- unique_ids(fam[[2L]], fams[f])
    rows <- match(sample_ids, physical_ids)
    if (anyNA(match(available, physical_ids))) .gsim_stop("Glist sample IDs absent from FAM: ", fams[f])
    bim <- read_metadata(bims[f]); physical_markers <- unique_ids(bim[[2L]], bims[f])
    mapping <- match(as.character(Glist$rsids[[f]]), physical_markers)
    if (anyNA(mapping)) .gsim_stop("Glist marker IDs absent from BIM: ", bims[f])
    if (!is.null(Glist$mchr) && (length(Glist$mchr) != length(beds) ||
        is.na(Glist$mchr[f]) || Glist$mchr[f] != length(mapping)))
      .gsim_stop("Glist$mchr must match the per-file metadata catalog.")
    for (allele in c("a1", "a2")) {
      if (!is.null(Glist[[allele]])) {
        expected <- Glist[[allele]][[f]]
        observed <- bim[[if (allele == "a1") 5L else 6L]][mapping]
        if (!identical(as.character(expected), observed))
          .gsim_stop("Glist$", allele, " disagrees with physical BIM orientation: ", bims[f])
      }
    }
    local <- match(marker_ids, as.character(Glist$rsids[[f]]))
    hit <- which(!is.na(local))
    file_index[hit] <- f
    variant_index[hit] <- mapping[local[hit]]
    files[[f]] <- list(enc2utf8(beds[f]), as.integer(nrow(fam)), as.integer(nrow(bim)), as.integer(rows))
    # Validate every declared BED's extent/header before sampling effects. This
    # bounded read is metadata validation, not a whole-panel scan.
    invisible(.Call(C_gsim_bed_read_selected, files[[f]][[1L]], files[[f]][[2L]],
                    files[[f]][[3L]], 1L, 1L))
  }
  list(files = files, file_index = file_index, variant_index = variant_index,
       marker_ids = marker_ids, sample_ids = sample_ids)
}

.gsim_bed_getG <- function(plan) {
  function(Glist, rsids, ids, chr = NULL, impute = TRUE, scale = FALSE) {
    pos <- match(rsids, plan$marker_ids)
    if (anyNA(pos) || !identical(as.character(ids), plan$sample_ids) || scale)
      .gsim_stop("BED block does not match the validated selection.")
    W <- matrix(NA_real_, length(ids), length(pos), dimnames = list(ids, rsids))
    for (f in unique(plan$file_index[pos])) {
      take <- which(plan$file_index[pos] == f)
      spec <- plan$files[[f]]
      W[, take] <- .Call(C_gsim_bed_read_selected, spec[[1L]], spec[[2L]], spec[[3L]],
                         spec[[4L]], plan$variant_index[pos[take]])
    }
    W
  }
}

.gsim_glist_statistics <- function(Glist, marker_ids, sample_ids, chr_map,
                                    getG_fun, chunk_size, standardize) {
  means <- replacements <- sds <- maf <- numeric(length(marker_ids))
  valid_dosage <- TRUE
  # Preserve canonical causal order; group only within each bounded read for
  # callbacks that require a chromosome argument.
  for (start in seq.int(1L, length(marker_ids), by = chunk_size)) {
    idx <- start:min(start + chunk_size - 1L, length(marker_ids))
    raw <- .gsim_load_glist_markers(Glist, marker_ids[idx], sample_ids,
                                   chr_map, getG_fun, chunk_size)
    valid_dosage <- valid_dosage && all(raw[is.finite(raw)] >= 0 & raw[is.finite(raw)] <= 2)
    maf[idx] <- .gsim_maf_from_W(raw)
    replacements[idx] <- colMeans(raw, na.rm = TRUE)
    W <- .gsim_impute_and_standardize(raw, standardize = FALSE)
    means[idx] <- colMeans(W)
    sds[idx] <- if (standardize) apply(W, 2L, stats::sd) else 1
    if (any(!is.finite(sds[idx]) | sds[idx] <= 0))
      .gsim_stop("Causal genotype columns have zero/non-finite variance: ",
                 paste(marker_ids[idx][!is.finite(sds[idx]) | sds[idx] <= 0], collapse = ", "))
  }
  if (!valid_dosage) maf[] <- NA_real_
  list(maf = maf, replacement = replacements, center = if (standardize) means else means * 0,
       scale = sds)
}

.gsim_glist_accumulate <- function(plan, statistics, B, Glist, marker_ids,
                                    sample_ids, chr_map, getG_fun, chunk_size,
                                    standardize) {
  if (!is.null(plan)) {
    pos <- match(marker_ids, plan$marker_ids)
    return(.Call(C_gsim_bed_accumulate, plan$files, plan$file_index[pos],
                 plan$variant_index[pos], statistics$replacement,
                 statistics$center, statistics$scale, B))
  }
  # Explicit custom storage adapters must be deterministic and honor bounded
  # requests. No implicit qgg/dense fallback exists for unsupported storage.
  G <- matrix(0, length(sample_ids), ncol(B))
  for (start in seq.int(1L, length(marker_ids), by = chunk_size)) {
    idx <- start:min(start + chunk_size - 1L, length(marker_ids))
    W <- .gsim_load_glist_markers(Glist, marker_ids[idx], sample_ids,
                                 chr_map, getG_fun, chunk_size)
    W <- .gsim_impute_and_standardize(W, standardize)
    G <- G + W %*% B[idx, , drop = FALSE]
  }
  G
}
