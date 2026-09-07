#' Generate gene sets conditional on known causal SNPs
#'
#' Construct controlled gene-set compositions for prioritization experiments.
#' A gene is causal when at least one supplied causal SNP maps to it. This is
#' structural truth, not a significance or enrichment label. No SNP effects,
#' phenotypes, or enrichment statistics are generated or modified.
#'
#' @param snp_gene_map Data frame with character `snp` and `gene` columns.
#'   IDs must be nonmissing and nonempty and are matched exactly. Duplicate
#'   pairs are removed; a SNP may map to multiple genes. Represent intergenic
#'   SNPs by absent mappings, not the reserved gene label `"intergenic"`.
#'   All unique mapped genes form the universe; genes with no mapped SNPs
#'   cannot be represented by this interface. Extra columns are ignored.
#' @param causal_snps Unique character vector of known causal SNP IDs, normally
#'   `simulation$causal_rsids`. `character(0)` permits all-null experiments.
#'   IDs absent from the map are returned and produce a concise warning.
#' @param set_size Positive integer counts of genes per set.
#' @param n_causal Nonnegative integer counts of causal genes per set, no larger
#'   than `set_size`. Equal-length vectors define paired scenarios; either
#'   scalar may be broadcast. Other length mismatches and infeasible pool sizes
#'   are errors. There is no implicit Cartesian product.
#' @param n_sets One positive integer: replicates per scenario.
#' @param seed Optional seed, passed to `set.seed()` at entry as in `gsim()`.
#'   The global RNG state is advanced, not restored. NULL uses the current RNG
#'   stream. Reproducibility assumes the same R RNG kind and version.
#'
#' @details Genes are sorted in UTF-8 radix order before sampling. Each set
#'   samples uniformly without replacement from the causal pool, then from the
#'   noncausal pool; its combined members are returned in canonical order.
#'   Pools are reused across sets, so natural overlap and identical replicates
#'   are allowed. Mapping-row and causal-SNP order do not affect seeded results.
#'   Set IDs are `scenario_<i>_replicate_<j>`, in scenario-major order.
#'
#'   Storage scales with mappings, genes, scenarios, and generated memberships;
#'   no dense gene-by-set or SNP-by-set matrix is allocated. If deriving SNP
#'   memberships downstream, deduplicate shared SNPs within each set. Multiple
#'   annotations must never multiply a SNP's supplied effect.
#'
#' @return A serializable ordinary list:
#' \describe{
#'   \item{sets}{Named list of unique gene-ID vectors.}
#'   \item{membership}{Data frame with `set_id`, `gene`, and logical `causal`.}
#'   \item{truth}{One row per set: `set_id`, integer `scenario` and `replicate`,
#'     `requested_size`, `requested_n_causal`, `realized_size`,
#'     `realized_n_causal`, and `causal_fraction`.}
#'   \item{genes}{Canonical `gene`, logical `causal`, and integer `n_snps` and
#'     `n_causal_snps`, counting unique SNPs per gene.}
#'   \item{unmapped_causal_snps}{Canonical causal SNP IDs absent from the map.}
#'   \item{settings}{Seed, paired `scenarios`, `n_sets`, and ordering, sampling,
#'     overlap and RNG semantics.}
#' }
#' @references Gholipourshahraki et al. (2024). Evaluation of Bayesian Linear
#'   Regression models for gene set prioritization in complex diseases.
#'   PLoS Genetics 20(11): e1011463. \doi{10.1371/journal.pgen.1011463}.
#' @examples
#' mapping <- data.frame(snp = paste0("s", 1:6), gene = paste0("g", 1:6))
#' sets <- gsim_gene_sets(mapping, c("s1", "s2"), set_size = 3,
#'                        n_causal = 0:2, n_sets = 2, seed = 41)
#' sets$truth
#' @export
gsim_gene_sets <- function(snp_gene_map, causal_snps, set_size, n_causal,
                           n_sets = 1L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  check_ids <- function(x, label) {
    if (!is.character(x) || !is.null(dim(x)) || anyNA(x) || any(!nzchar(x)))
      stop(label, " must contain nonmissing, nonempty character IDs.", call. = FALSE)
    unname(enc2utf8(x))
  }
  counts <- function(x, label, minimum) {
    if (!is.numeric(x) || !is.null(dim(x)) || !length(x) ||
        any(!is.finite(x)) || any(x != floor(x)) ||
        any(x < minimum | x > .Machine$integer.max))
      stop(label, " must contain finite integer counts between ", minimum,
           " and .Machine$integer.max.", call. = FALSE)
    unname(as.integer(x))
  }
  if (!is.data.frame(snp_gene_map) ||
      sum(names(snp_gene_map) == "snp") != 1L ||
      sum(names(snp_gene_map) == "gene") != 1L)
    stop("snp_gene_map must be a data.frame with snp and gene columns.", call. = FALSE)
  snp <- check_ids(snp_gene_map$snp, "snp_gene_map$snp")
  gene <- check_ids(snp_gene_map$gene, "snp_gene_map$gene")
  if (any(gene == "intergenic"))
    stop("Represent intergenic SNPs by absent mappings, not gene 'intergenic'.", call. = FALSE)
  causal_snps <- check_ids(causal_snps, "causal_snps")
  if (anyDuplicated(causal_snps)) stop("causal_snps must be unique.", call. = FALSE)
  causal_snps <- sort(causal_snps, method = "radix")
  set_size <- counts(set_size, "set_size", 1L)
  n_causal <- counts(n_causal, "n_causal", 0L)
  n_sets <- counts(n_sets, "n_sets", 1L)
  if (length(n_sets) != 1L) stop("n_sets must be one positive integer.", call. = FALSE)
  scenarios_n <- max(length(set_size), length(n_causal))
  if (!(length(set_size) %in% c(1L, scenarios_n)) ||
      !(length(n_causal) %in% c(1L, scenarios_n)))
    stop("set_size and n_causal must have equal lengths or one must be scalar.", call. = FALSE)
  set_size <- rep(set_size, length.out = scenarios_n)
  n_causal <- rep(n_causal, length.out = scenarios_n)
  mapping <- unique(data.frame(snp = snp, gene = gene, stringsAsFactors = FALSE))
  universe <- sort(unique(gene), method = "radix")
  gene_index <- match(mapping$gene, universe)
  causal_mapping <- mapping$snp %in% causal_snps
  causal_counts <- tabulate(gene_index[causal_mapping], nbins = length(universe))
  genes <- data.frame(gene = universe, causal = causal_counts > 0L,
                      n_snps = tabulate(gene_index, nbins = length(universe)),
                      n_causal_snps = causal_counts, stringsAsFactors = FALSE)
  causal_pool <- universe[genes$causal]
  null_pool <- universe[!genes$causal]
  bad <- which(n_causal > set_size | n_causal > length(causal_pool) |
                 set_size - n_causal > length(null_pool))
  if (length(bad)) {
    i <- bad[1L]
    stop("Scenario ", i, " requests size ", set_size[i], " and ", n_causal[i],
         " causal genes; available: ", length(causal_pool), " causal and ",
         length(null_pool), " noncausal genes. Require 0 <= n_causal <= set_size and sufficient pools.",
         call. = FALSE)
  }
  if (sum(as.double(set_size)) * n_sets > .Machine$integer.max)
    stop("Requested memberships exceed the supported data.frame row limit.", call. = FALSE)
  unmapped <- setdiff(causal_snps, mapping$snp)
  if (length(unmapped)) warning(length(unmapped),
    " causal SNP(s) absent from snp_gene_map; see unmapped_causal_snps.", call. = FALSE)
  scenarios <- data.frame(scenario = seq_len(scenarios_n), set_size = set_size,
                          n_causal = n_causal)
  scenario <- rep(seq_len(scenarios_n), each = n_sets)
  replicate <- rep(seq_len(n_sets), times = scenarios_n)
  set_ids <- paste0("scenario_", scenario, "_replicate_", replicate)
  draw <- function(pool, size) {
    if (!size) return(character(0))
    pool[sample.int(length(pool), size, replace = FALSE)]
  }
  sets <- lapply(scenario, function(i) {
    sort(c(draw(causal_pool, n_causal[i]),
           draw(null_pool, set_size[i] - n_causal[i])), method = "radix")
  })
  names(sets) <- set_ids
  members <- unlist(sets, use.names = FALSE)
  flags <- genes$causal[match(members, universe)]
  membership <- data.frame(set_id = rep(set_ids, lengths(sets)), gene = members,
                            causal = flags, stringsAsFactors = FALSE)
  realized <- vapply(sets, function(x) sum(x %in% causal_pool), integer(1))
  truth <- data.frame(set_id = set_ids, scenario = scenario, replicate = replicate,
    requested_size = set_size[scenario], requested_n_causal = n_causal[scenario],
    realized_size = unname(lengths(sets)), realized_n_causal = unname(realized),
    causal_fraction = unname(realized / lengths(sets)), stringsAsFactors = FALSE)
  list(sets = sets, membership = membership, truth = truth, genes = genes,
    unmapped_causal_snps = unmapped,
    settings = list(seed = seed, scenarios = scenarios, n_sets = n_sets,
      ordering = "UTF-8 radix genes; scenario-major sets",
      sampling = "uniform without replacement within each causal/noncausal pool per set",
      overlap = "pools reused across sets; overlap and identical sets allowed",
      rng = "set.seed at entry when supplied; global state advanced, not restored"))
}
