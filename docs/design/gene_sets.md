# Controlled gene-set simulation

## Scope and scientific interpretation

`gsim_gene_sets()` constructs sets conditional on supplied causal SNP truth.
The design follows the gene-set simulation experiment in Gholipourshahraki
et al. (2024), [Evaluation of Bayesian Linear Regression models for gene set
prioritization in complex diseases](https://doi.org/10.1371/journal.pgen.1011463),
PLoS Genetics 20(11): e1011463, specifically its “Simulation of gene sets” section.
The paper motivates varying total and causal gene counts, replicates and
zero-causal controls. This implementation provides that set-construction design;
it does not reproduce the paper's full simulation or evaluate its methods.

A gene is causal if at least one supplied causal SNP maps to it. This is structural
simulation truth, not statistical significance. Nonzero causal content does not
by itself label a set statistically enriched: that depends on background and
analysis. Effects remain in the original `gsim()` result and are never redrawn,
rescaled or multiplied by annotation membership here.

## Input and sampling contract

```r
gsim_gene_sets(snp_gene_map, causal_snps, set_size, n_causal,
               n_sets = 1L, seed = NULL)
```

The map is a data.frame with character `snp` and `gene` columns; extra columns
are ignored. IDs must be nonmissing and nonempty and match exactly, without
trimming or case folding. Repeated identical pairs are deduplicated. One SNP
may map to several genes, making each causal if that SNP is causal. Per-gene
SNP counts count distinct SNPs, never repeated rows. `causal_snps` must be unique;
`character(0)` supports all-null experiments.

All unique mapped genes form the universe. Genes with no mapped SNPs are outside
this interface. Intergenic SNPs have absent mappings; the fictitious gene label
`intergenic` is rejected. Unknown causal SNPs are returned in
`unmapped_causal_snps` with a concise count warning. They do not create genes.

`set_size` and `n_causal` are paired finite integer vectors. Equal lengths define
one scenario per pair; one scalar may broadcast to the other vector. Other
length mismatches fail; no Cartesian product is formed. `n_sets` is one positive
integer specifying replicates per scenario. Every size is at least one and
0 <= n_causal <= set_size. Insufficient causal or noncausal pools fail with the
scenario number, requested counts and both available pool counts. Counts are
limited to R's integer range, and total memberships to its ordinary integer
row-count range. Requests are never capped or altered to fit the pools.

For each set, sample uniformly without replacement from the causal pool and then
from the noncausal pool, taking exactly the requested number from each. Combined
members have no duplicates. Pools are reused unchanged across all sets, allowing
natural overlap and identical replicate sets. No overlap rate is specified.
If downstream SNP memberships are derived, take unique SNP IDs within each set:
shared SNPs and multiple annotations must never multiply their supplied effects.

## Determinism and returned truth

Genes use canonical UTF-8 radix ordering, independent of mapping row order;
causal SNPs are canonicalized too. Returned gene vectors use that same order.
Set IDs are `scenario_<i>_replicate_<j>`, in scenario-major order. Reordering
mapping rows or causal SNP IDs preserves the complete result for the same seed.
Scenario order itself is meaningful and is not canonicalized.

As in `gsim()`, a supplied seed calls `set.seed()` at entry, before validation,
and subsequent draws advance the global R RNG state; the old state is not
restored. NULL continues the current stream. No RNG kind is changed. Exact
reproducibility assumes the same R version and RNG kind. Zero-count pool draws
consume no random values; nonzero draws use `sample.int()` without replacement.

The ordinary serializable list contains:

- `sets`: named gene-ID vectors.
- `membership`: one row per set–gene pair, with `set_id`, `gene`, and `causal`.
- `truth`: set ID, scenario, replicate, requested/realized size and causal count,
  and realized causal fraction.
- `genes`: canonical gene IDs, causal flags, unique mapped SNP counts (`n_snps`)
  and unique causal SNP counts (`n_causal_snps`).
- `unmapped_causal_snps`: sorted causal SNP IDs absent from the map.
- `settings`: seed, paired scenarios, replicate count and sampling/overlap/RNG
  conventions.

The [public function reference](../../man/gsim_gene_sets.Rd) defines column names.
Storage is proportional to input mapping pairs, genes, scenarios, sets and
emitted memberships, including temporary deduplication, indexing and sampling
vectors. Gene vectors and the long membership table both retain memberships;
there is no dense gene × set or SNP × set matrix. This is not constant-memory
output, and large requested membership counts still require commensurate memory.

## Example and focused verification

Run the [complete local example](../../inst/examples/gene_sets.R) after installing
gsim:

```r
source("inst/examples/gene_sets.R")
```

It simulates 64 individuals and 24 markers with four causal SNPs, constructs an
explicitly synthetic one-SNP/one-gene map, then creates size-12 sets with 0, 2
and 4 causal genes and three replicates. Two size-12 controls from twenty
noncausal genes necessarily overlap. The script prints truth and shared genes.
Substitute a real two-column SNP–gene map with SNP IDs matching
`simulation$causal_rsids`; no coordinate generation or annotation download is
performed. Existing `gsim()` alone remains responsible for effects/phenotypes.

[Focused tests](../../tests/testthat/test-gene-sets.R) cover exact composition,
null controls, empty causal input, many-gene mappings, pair deduplication,
infeasible inputs, unmapped reporting, canonical ordering, independent reference
RNG draws, within-set uniqueness, reused pools, serialization and actual gsim
causal output. A compact 24-gene storage fixture compares 20 and 40 sets of four
genes (80 and 160 membership rows), without constructing an incidence matrix.
These checks are set-generation qualification, not enrichment performance.

Local verification on 2026-09-07 used R 4.4.1 on Windows and one isolated 0.16.0
installation. All 108 focused checks passed (no test failures, warnings or skips),
and the complete example ran successfully. A 24-gene fixture with SNP IDs s1
through s24, genes g1 through g24, causal SNPs s1 through s4, seed 9, size 4 and
one causal gene occupied 24,392 bytes for 20 sets / 80 memberships and 37,816
bytes for 40 sets / 160 memberships (`object.size`, including names and R
overhead; not peak process memory). Entry was clean `main` at
`3f540c6f70862ea6bd87617658a44350edd71584`, version 0.15.0. Existing phenotype
and native simulation code were unchanged; no full suite or benchmark ran.

Deferred capabilities: synthetic gene placement, controlled overlap rates,
real-annotation acquisition, and enrichment-method evaluation. No enrichment
testing, Bayesian inference or additional marker-effect model is included.
