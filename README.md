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

## Provenance

The initial implementation was extracted from `sblr` commit
`e9532f8b852f973f34f531a1cc9101da75e1f0ad`, using `R/gsim.R`,
`R/gsim_internal.R`, and `tests/testthat/test-gsim.R` as the canonical source
paths.
