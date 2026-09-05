# gsim

`gsim` is a compact R package for simulating genomic phenotypes with exact
marker-level truth for validation and methodological studies. Genotypes can be
simulated internally as independent binomial markers, supplied as an in-memory
matrix, or read from an optional `qgg::Glist` when `qgg` is installed.

The simulator supports BayesC, BayesR, major-plus-polygenic, MAF-dependent,
clustered, and fixed-effect architectures. It can simulate one or multiple
traits, annotation-informed component probabilities, and marker-specific
active-effect variance multipliers. Optional marginal summary statistics can be
generated from the same simulated phenotype.

Returned objects include phenotype, genetic and residual components, exact
marker effects and states, causal markers, probability surfaces, targets,
realized quantities, and settings/provenance used to generate the data. These
outputs are intended to make validation studies directly auditable.

## Examples

Simulate genotypes and a phenotype:

```r
library(gsim)

sim <- gsim(n = 200, m = 500, n_causal = 20, seed = 1)
sim
```

Simulate from a caller-provided genotype matrix:

```r
W <- matrix(rbinom(150 * 100, size = 2, prob = 0.3), 150, 100)
colnames(W) <- paste0("m", seq_len(ncol(W)))
sim <- gsim(W = W, architecture = "bayesr", n_causal = 10, seed = 2)
```

Real `qgg::Glist` inputs are supported when `qgg` is installed. Scientific
validation studies remain in packages such as `sblrbench`; they are not bundled
with `gsim`.

## Packed reference workflow

The supported packed workflow imports a strict phased biallelic text VCF and
then simulates HAPNEST-compatible founders plus Mendelian pedigree descendants
without a dense whole-genome allele or genotype matrix:

```r
reference <- gsim_import_vcf(
  "reference.vcf", genetic_map, "reference"
)
result <- gsim_simulate(
  reference, pedigree, populations, ancestry_weights, mutation_age,
  N, Ne, rho, seed = 123, output = "simulation", format = "hap"
)
```

The optional native `gbits` and `gmat` libraries are located through
`GSIM_GBITS_LIBRARY` and `GSIM_GMAT_LIBRARY`. VCF REF is bit 0/BIM A2, ALT is
bit 1/BIM A1, and phased GT left/right order is retained as H1/H2. The importer
accepts only uncompressed VCF, complete phased diploid GT, and uppercase
biallelic A/C/G/T SNPs. The required map supplies exact cumulative cM values;
there is no interpolation. Each chromosome is loaded, simulated, written, and
released before the next.

## Pedigree and record workloads

`gsim_pedigree()` creates scalable multigenerational pedigree domains with
restricted sire and dam pools, overlapping generations, missing parents, later
founders, deliberately unphenotyped animals, and separate canonical and external
orders. Parent-before-offspring ordering is explicit in `canonical_order`; the
returned pedigree table uses the reproducibly arbitrary `external_order`.

`gsim_pedigree_records()` turns one pedigree into one selected model view:
single-trait, two-trait with incomplete observation patterns, or irregular
longitudinal random regression. Every observed phenotype is one scalar row.
Missing traits and times are absent records, not imputed values. The fixed design
is a sorted one-based triplet list (`row`, `column`, `value`) with four stored
entries per observation, so no large dense incidence matrix is returned.

Longitudinal views store a basis row aligned with every observed record. Optional
prediction records provide animal, new time, basis, and fitted truth without an
observed residual or phenotype.

The pedigree latent values are deterministic solver-workload values formed by a
scalable parent-average recursion. They are not an inbreeding-aware exact draw
from a numerator-relationship covariance, and covariance recovery is therefore
not a validation target. Their purpose is to supply identical model inputs and
right-hand sides for sparse solver parity studies.

```r
ped <- gsim_pedigree(
  n_generations = 5, animals_per_generation = 40,
  sires_per_generation = 6, dams_per_generation = 12, seed = 10
)
long <- gsim_pedigree_records(
  ped, model = "longitudinal", prediction_records = TRUE, seed = 11
)
```

The non-test script `tools/qualification/pedigree_solver_workload.R` contains the
fixed 50,000-animal construction qualification. It prints counts, dimensions,
object sizes, elapsed construction times, and deterministic checksums, produces
no permanent output by default, and is never run by `R CMD check`.

## Provenance

The initial implementation was extracted from `sblr` commit
`e9532f8b852f973f34f531a1cc9101da75e1f0ad`, using `R/gsim.R`,
`R/gsim_internal.R`, and `tests/testthat/test-gsim.R` as the canonical source
paths.
