# gsim

`gsim` simulates genotypes, pedigrees, and genomic phenotypes for validation
and methodological studies. Its standalone packed workflow imports phased VCF
reference panels or opens HAP/BIM/FAM, samples unrelated synthetic founders
using a HAPNEST-informed model, and transmits their chromosomes through a
pedigree by Mendelian meiosis. Reusable phased HAP and unphased BED outputs
connect genotype simulation to phenotype studies with exact marker-level truth.
Phenotype inputs can also be independent binomial markers, a caller-provided
matrix, or a `qgg::Glist`. qgg is optional generally and required for the
currently supported Glist workflow.

See the [documentation index](docs/README.md) for scientific contracts, public
API references, complete examples, and qualification evidence.

The simulator supports BayesC, BayesR, major-plus-polygenic, MAF-dependent,
clustered, and fixed-effect architectures. It can simulate one or multiple
traits, annotation-informed component probabilities, direct marker-specific
causal probabilities, and marker-specific active-effect variance weights.
Optional marginal summary statistics can be generated from the same simulated
phenotype.

Returned objects include phenotype, genetic and residual components, exact
marker effects and states, causal markers, probability surfaces, targets,
realized quantities, and settings/provenance used to generate the data. These
outputs are intended to make validation studies directly auditable.

## Installation

Install the current release directly from GitHub with `remotes`:

```r
remotes::install_github("psoerensen/gsim")
```

`remotes` is used only as an installation helper and is not a `gsim` package
dependency. From a local checkout, the equivalent base-R command is:

```text
R CMD INSTALL --preclean .
```

For RStudio development, open `gsim.Rproj` and use **Build > Install** (or the
Install button in the Build pane). The project is configured as an R package
project and invokes base `R CMD INSTALL --preclean`; it does not require
`devtools`, add RStudio as a package dependency, or include RStudio's
`.Rproj.user` state in Git or built source packages.

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

Membership and active-effect variance can be specified independently. A direct
probability is a Bernoulli non-null probability; a multiplier is a relative
active-effect variance weight:

```r
marker_ids <- colnames(W)
q <- setNames(seq(0.02, 0.20, length.out = length(marker_ids)), marker_ids)
w <- setNames(exp(seq(-0.5, 0.5, length.out = length(marker_ids))), marker_ids)
sim <- gsim(
  W = W, architecture = "bayesr",
  causal_probability = q, marker_multipliers = w, seed = 3
)
```

Real `qgg::Glist` inputs are supported when `qgg` is installed. Scientific
validation studies remain in packages such as `sblrbench`; they are not bundled
with `gsim`.

When a Glist contains `maf` and `ldscores` metadata, conditional effect-variance
weights can be derived without loading the complete genotype matrix. The scalar
annotation score is separate from the SBayesRC annotation matrix `A`.
This workflow sketch requires a prepared `Glist` with aligned metadata:

```r
marker_ids <- as.character(unlist(Glist$rsidsLD, use.names = FALSE))
q <- setNames(rep(0.05, length(marker_ids)), marker_ids)
s <- setNames(rep(1, length(marker_ids)), marker_ids)
sim <- gsim(
  Glist = Glist, architecture = "bayesr",
  causal_probability = q,
  a = -0.4, b = -1, c = 0.5, annotation_score = s,
  seed = 4
)
```

This illustrates a HAPNEST-like variance parameterization, not a separate
architecture or a recommended universal default. `causal_probability` controls
non-null membership; active mixture proportions remain conditional on `pi`;
`a`, `b`, `c`, and `annotation_score` control variance only after a marker is
causal. Effects scale by the square root of the resulting variance weight.

## Packed reference workflow

The supported packed workflow explicitly separates a real reference panel, an
unrelated synthetic base population, and its Mendelian pedigree descendants.
The packed genotype stages never construct a dense whole-genome allele or
genotype matrix. This workflow sketch requires a VCF, genetic map, named model
inputs, and a pedigree whose founder IDs match the generated base IDs:

```r
reference <- gsim_import_vcf("reference.vcf.gz", genetic_map, "reference")
base <- gsim_simulate_founders(
  reference, n = 10000, populations = populations,
  ancestry_weights = ancestry_weights, mutation_age = mutation_age,
  N = N, Ne = Ne, rho = rho, seed = 123, output = "base",
  batch_size = 4096, threads = 8
)
result <- gsim_simulate_pedigree(
  base, pedigree, seed = 456, output = "pedigree", format = "hap"
)
```

All packed storage, VCF parsing, and dataset writing code is compiled directly
into `gsim`; no sibling packages, shared-library paths, or environment variables
are required. VCF REF is bit 0/BIM A2, ALT is bit 1/BIM A1, and phased GT
left/right order is retained as H1/H2. The importer
retains complete phased diploid GT at uppercase biallelic A/C/G/T SNPs.
Unsupported biological records may be counted and skipped. A sparse physical
map supplies cumulative cM knots, with deterministic interpolation and no
extrapolation. Founder batches are aligned to 64-sample packed words and written
directly into final marker-major HAP positions. Static native workers own
disjoint output words, so batch size and thread scheduling do not change output.
The founder seed affects only the base population; the pedigree seed affects
only meiosis. HAP retains phase and is reusable across pedigree runs, while BED
is an unphased dosage output. Pedigree meiosis remains single-threaded because
parents and children can occupy the same mutable packed word.

See [the direct 1000 Genomes chromosome 22 example](inst/examples/1000G_chr22.R)
for an internet-enabled GRCh37 workflow using the official approximately
196 MB IGSR VCF directly, without bcftools or htslib.

## Genotypes to phenotypes

Packed simulation output can enter the existing phenotype engine without a new
adapter. This workflow sketch requires prepared reference and pedigree inputs;
`...` stands for required founder-model arguments, not runnable R code:

```r
reference <- gsim_reference("reference")
base <- gsim_simulate_founders(reference, ..., output = "base")
pedigree_bed <- gsim_simulate_pedigree(
  base, pedigree, seed = 2, output = "pedigree", format = "bed"
)
Glist <- qgg::gprep(
  study = "simulation",
  bedfiles = pedigree_bed$paths[["bed"]],
  bimfiles = pedigree_bed$paths[["bim"]],
  famfiles = pedigree_bed$paths[["fam"]]
)
phenotype <- gsim(Glist = Glist, n_causal = 20, seed = 3)
```

BED dosage is H1 + H2, BIM order becomes `Glist$rsids`, FAM order becomes
`Glist$ids`, and `gprep()` supplies `Glist$maf`. LD scores must be prepared by
the established Glist workflow or supplied as a complete named vector. During
phenotype simulation, gsim selects causal markers first and requests only those
columns from `qgg::getG()` unless summary statistics are requested. Selected
causal columns are decoded into a dense R matrix; this is not a chromosome-local
packed phenotype engine. The supported Glist route uses qgg; compatibility with
gsuite-generated Glist objects is not established by the retained evidence.

The three marker-level controls are distinct: `q_j` is causal probability,
active `pi_k` values are conditional mixture proportions, and `w_j` is the
conditional effect-variance multiplier. See the
[complete local-data example](inst/examples/end_to_end_phenotype.R)
(run after `library(gsim)`, with qgg installed) and the
[qualification contract](docs/qualification/genotypes_to_phenotypes.md) for the
default, causal-probability, variance-weight, and combined configurations.

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

HAPNEST scientific attribution and gbits/gmat component provenance are retained
in the [packed workflow guide](docs/design/packed_chromosome_simulation.md) and
the applicable [copyright and license notices](inst/COPYRIGHTS).
