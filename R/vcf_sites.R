#' Stream VCF site metadata in batches
#'
#' Reads plain, gzip or BGZF VCF sequentially without decoding sample
#' genotypes. BCF and indexed access are unsupported.
#'
#' @param vcf Path to a VCF file.
#' @param callback Function called with each data frame of site records. When
#'   writing a subset, it must return one nonmissing logical per input row.
#' @param info Character vector of INFO keys. Values are returned as text,
#'   absent keys as NA, and flags as "TRUE".
#' @param chromosome Optional exact chromosome label.
#' @param region Optional inclusive positive integer bounds, requiring
#'   a chromosome.
#' @param batch_size Maximum number of records passed to the callback at once.
#' @param output Optional path for a new gzip VCF containing the original header
#'   and selected records.
#' @return Invisibly, the VCF header lines.
#' @details Callback data frames contain numeric `pos` and character
#'   `chrom`, `id`, `ref`, `alt`, `qual`, and `filter` columns, followed
#'   by the requested INFO values in columns prefixed with `info_`.
#' @export
gsim_scan_vcf <- function(vcf, callback, info = character(), chromosome = NULL,
                          region = NULL, batch_size = 4096L, output = NULL) {
  if (!is.character(vcf) || length(vcf) != 1L || is.na(vcf) || !nzchar(vcf)) {
    stop("vcf must be one path", call. = FALSE)
  }
  if (!is.function(callback) || !is.character(info) || anyNA(info) ||
      any(!nzchar(info)) || anyDuplicated(info)) {
    stop("Supply a callback and unique nonempty INFO keys", call. = FALSE)
  }
  if (!is.null(chromosome) &&
      (!is.character(chromosome) || length(chromosome) != 1L ||
       is.na(chromosome) || !nzchar(chromosome))) {
    stop("chromosome must be one nonempty label", call. = FALSE)
  }
  if (!is.null(region) &&
      (is.null(chromosome) || !is.numeric(region) || length(region) != 2L ||
       any(!is.finite(region) | region < 1 | region != floor(region)) ||
       region[1L] > region[2L])) {
    stop("region must be positive integer bounds with a chromosome", call. = FALSE)
  }
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      !is.finite(batch_size) || batch_size < 1 || batch_size > .Machine$integer.max ||
      batch_size != floor(batch_size)) {
    stop("batch_size must be a positive integer", call. = FALSE)
  }

  reader <- .Call(C_gsim_vcf_sites_open, enc2utf8(normalizePath(vcf, mustWork = TRUE)))
  pointer <- reader[[1L]]
  header <- reader[[2L]]
  on.exit(.Call(C_gsim_vcf_sites_close, pointer), add = TRUE)
  if (!is.null(output)) {
    if (!is.character(output) || length(output) != 1L || is.na(output) ||
        !nzchar(output) || file.exists(output)) {
      stop("output must be a new file path", call. = FALSE)
    }
    lock <- paste0(output, ".lock")
    if (!dir.create(lock, showWarnings = FALSE)) {
      stop("Could not create a staging directory for the VCF subset", call. = FALSE)
    }
    on.exit(unlink(lock, recursive = TRUE, expand = FALSE), add = TRUE)
    if (file.exists(output)) stop("output must be a new file path", call. = FALSE)
    staged <- file.path(lock, "subset.vcf.gz")
    connection <- gzfile(staged, open = "wb", compression = 1L)
    on.exit(if (!is.null(connection)) close(connection), add = TRUE, after = FALSE)
    writeLines(header, connection, useBytes = TRUE)
  }
  columns <- c("pos", "chrom", "id", "ref", "alt", "qual", "filter",
               if (length(info)) paste0("info_", info))
  repeat {
    batch <- .Call(C_gsim_vcf_sites_next, pointer, as.integer(batch_size),
                   enc2utf8(info), !is.null(output),
                   if (is.null(chromosome)) "" else enc2utf8(chromosome),
                   as.double(region))
    n <- batch[[2L]]
    if (!n) break
    rows <- seq_len(n)
    records <- lapply(batch[[1L]][seq_along(columns)], function(value) value[rows])
    names(records) <- columns
    keep <- callback(as.data.frame(records, stringsAsFactors = FALSE, optional = TRUE))
    if (!is.null(output)) {
      if (!is.logical(keep) || length(keep) != n || anyNA(keep)) {
        stop("callback must return one nonmissing logical per input row", call. = FALSE)
      }
      writeLines(batch[[1L]][[length(columns) + 1L]][rows][keep], connection, useBytes = TRUE)
    }
  }
  if (!is.null(output)) {
    close(connection)
    connection <- NULL
    if (file.exists(output) || !file.rename(staged, output)) {
      stop("Could not publish the VCF subset", call. = FALSE)
    }
  }
  invisible(header)
}
